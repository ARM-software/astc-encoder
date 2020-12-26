/* This pipeline is used for release testing, so it runs rarely.
 *
 * Test objectives for this pipeline are:
 *
 *   - Run the entire pipeline in less than 60 minutes.
 *   - Test builds on all supported operating systems.
 *   - Test builds on optimized compiler choices (i.e. prefer Clang over GCC).
 *   - Build only release variants.
 *   - Run full functional tests.
 *   - Run full image quality tests.
 *   - Code sign the binaries on supported operating systems.
 *   - Build the release package.
 *
 * The test matrix is not fully covered; e.g. we can assume compilers behave
 * similarly on different operating systems, so we test one compiler per OS.
 */

@Library('hive-infra-library@master') _

pipeline {
  agent none

  options {
    ansiColor('xterm')
    timestamps()
  }

  stages {
    stage('Build All') {
      parallel {
        /* Run static analysis on Linux */
        stage('Coverity') {
          agent {
            docker {
              image 'astcenc:2.4.0'
              registryUrl 'https://mobile-studio--docker.artifactory.geo.arm.com'
              registryCredentialsId 'cepe-artifactory-jenkins'
              label 'docker'
              alwaysPull true
            }
          }
          stages {
            stage('Clean') {
              steps {
                sh 'git clean -fdx'
              }
            }
            stage('Coverity') {
              steps {
                withCredentials([usernamePassword(credentialsId: 'jenkins-password',
                                                  usernameVariable: 'USERNAME',
                                                  passwordVariable: 'PASSWORD')]) {
                  sh script: '''#!/bin/bash
                    mkdir -p ${WORKSPACE}/occonfig

                    mkdir build_cov
                    cd build_cov

                    cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../ -DISA_AVX2=ON ..

                    cov-configure --template --compiler cc --comptype gcc
                    cov-configure --template --compiler c++ --comptype g++
                    cov-build --dir ${WORKSPACE}/intermediate make install
                    cov-analyze --dir ${WORKSPACE}/intermediate
                    cov-commit-defects --dir ${WORKSPACE}/intermediate \\
                                       --stream astcenc-master \\
                                       --url https://coverity.cambridge.arm.com \\
                                       --user jenkins@arm.com --password ${PASSWORD} \\
                                       --strip-path ${WORKSPACE}
                  '''
                }
              }
            }
          }
        }
        /* Build for Linux on x86-64 using Clang */
        stage('Linux') {
          agent {
            docker {
              image 'astcenc:2.3.0'
              registryUrl 'https://mobile-studio--docker.artifactory.geo.arm.com'
              registryCredentialsId 'cepe-artifactory-jenkins'
              label 'docker'
              alwaysPull true
            }
          }
          stages {
            stage('Clean') {
              steps {
                sh 'git clean -ffdx'
              }
            }
            stage('Build R') {
              steps {
                sh '''
                  export CXX=clang++-9
                  mkdir build_rel
                  cd build_rel
                  cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../ -DISA_AVX2=ON -DISA_SSE41=ON -DISA_SSE2=ON -DISA_NONE=ON ..
                  make install package -j4
                '''
              }
            }
            stage('Stash') {
              steps {
                dir('build_rel') {
                  stash name: 'astcenc-linux-x64', includes: '*.zip'
                  stash name: 'astcenc-linux-x64-hash', includes: '*.zip.sha256'
                }
              }
            }
            stage('Test') {
              steps {
                sh '''
                  python3 ./Test/astc_test_functional.py
                  python3 ./Test/astc_test_image.py --encoder=all --test-set Small
                '''
              }
            }
          }
        }
        /* Build for Windows on x86-64 using MSVC */
        stage('Windows') {
          agent {
            label 'Windows && x86_64'
          }
          stages {
            stage('Clean') {
              steps {
                bat 'git clean -ffdx'
              }
            }
            stage('Build R') {
              steps {
                bat '''
                  call c:\\progra~2\\micros~1\\2019\\buildtools\\vc\\auxiliary\\build\\vcvars64.bat
                  mkdir build_rel
                  cd build_rel
                  cmake -G "NMake Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../ -DISA_AVX2=ON -DISA_SSE41=ON -DISA_SSE2=ON ..
                  nmake install package
                '''
              }
            }
            stage('Stash') {
              steps {
                dir('build_rel') {
                  stash name: 'astcenc-windows-x64', includes: '*.zip'
                  stash name: 'astcenc-windows-x64-hash', includes: '*.zip.sha256'
                }
              }
            }
            stage('Test') {
              steps {
                bat '''
                  set Path=c:\\Python38;c:\\Python38\\Scripts;%Path%
                  call python ./Test/astc_test_image.py --test-set Small
                '''
              }
            }
          }
        }
        /* Build for macOS on x86-64 using Clang */
        stage('macOS') {
          agent {
            label 'mac && x86_64'
          }
          stages {
            stage('Clean') {
              steps {
                sh 'git clean -ffdx'
              }
            }
            stage('Build R') {
              steps {
                sh '''
                  mkdir build_rel
                  cd build_rel
                  cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../ -DISA_AVX2=ON -DISA_SSE41=ON -DISA_SSE2=ON ..
                  make install package -j1
                '''
              }
            }
            stage('Sign') {
              steps {
                dir('build_rel') {
                  withCredentials([sshUserPrivateKey(credentialsId: 'gerrit-jenkins-ssh',
                                                     keyFileVariable: 'SSH_AUTH_FILE')]) {
                    sh 'GIT_SSH_COMMAND="ssh -i $SSH_AUTH_FILE -o StrictHostKeyChecking=no" git clone ssh://eu-gerrit-1.euhpc.arm.com:29418/Hive/shared/signing'
                  }
                  withCredentials([usernamePassword(credentialsId: 'win-signing',
                                                    usernameVariable: 'USERNAME',
                                                    passwordVariable: 'PASSWORD')]) {
                    sh 'python3 ./signing/macos-client-wrapper.py ${USERNAME} *.zip'
                    sh 'rm -rf ./signing'
                  }
                }
              }
            }
            stage('Stash') {
              steps {
                dir('build_rel') {
                  stash name: 'astcenc-macos-x64', includes: '*.zip'
                  stash name: 'astcenc-macos-x64-hash', includes: '*.zip.sha256'
                }
              }
            }
            stage('Test') {
              steps {
                sh '''
                  export PATH=/usr/local/bin:$PATH
                  python3 ./Test/astc_test_image.py --test-set Small
                '''
              }
            }
          }
        }
      }
    }
    stage('Artifactory') {
      agent {
        docker {
          image 'astcenc:2.3.0'
          registryUrl 'https://mobile-studio--docker.artifactory.geo.arm.com'
          registryCredentialsId 'cepe-artifactory-jenkins'
          label 'docker'
          alwaysPull true
        }
      }
      options {
        skipDefaultCheckout true
      }
      stages {
        stage('Unstash') {
          steps {
            dir('upload/linux-x64') {
              unstash 'astcenc-linux-x64'
            }
            dir('upload/windows-x64') {
              unstash 'astcenc-windows-x64'
              withCredentials([sshUserPrivateKey(credentialsId: 'gerrit-jenkins-ssh',
                                                 keyFileVariable: 'SSH_AUTH_FILE')]) {
                sh 'GIT_SSH_COMMAND="ssh -i $SSH_AUTH_FILE -o StrictHostKeyChecking=no" git clone ssh://eu-gerrit-1.euhpc.arm.com:29418/Hive/shared/signing'
              }
              withCredentials([usernamePassword(credentialsId: 'win-signing',
                                                usernameVariable: 'USERNAME',
                                                passwordVariable: 'PASSWORD')]) {
                sh 'python3 ./signing/authenticode-client.py -u ${USERNAME} *.zip'
                sh 'rm -rf ./signing'
              }
            }
            dir('upload/macos-x64') {
              unstash 'astcenc-macos-x64'
            }
            dir('upload') {
              unstash 'astcenc-linux-x64-hash'
              // Don't keep Windows hash - we have signed binaries now
              // unstash 'astcenc-windows-x64-hash'
              unstash 'astcenc-macos-x64-hash'
            }
          }
        }
        stage('Upload') {
          steps {
            zip zipFile: 'astcenc.zip', dir: 'upload', archive: false
            cepeArtifactoryUpload(sourcePattern: 'astcenc.zip')
          }
        }
      }
      post {
        always {
          deleteDir()
        }
      }
    }
  }

  post {
    failure {
      script {
        slackSend channel: '#dsg-eng-astcenc', color: 'danger', message: "Build ${JOB_NAME} ${BUILD_NUMBER} failed. (<${BUILD_URL}|Open>)", teamDomain: 'arm-dsg', tokenCredentialId: 'jenkins-slack', username: 'jenkins'
      }
    }
  }
}
