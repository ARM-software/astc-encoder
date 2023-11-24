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

@Library('hive-infra-library@changes/86/295486/1') _

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
            kubernetes {
              yaml '''
apiVersion: v1
kind: Pod
spec:
  securityContext:
    runAsUser: 1000
    runAsGroup: 1000
  imagePullSecrets:
    - name: artifactory-ms-docker
  containers:
    - name: astcenc
      image: mobile-studio--docker.eu-west-1.artifactory.aws.arm.com/astcenc:3.2.0
      command:
        - sleep
      args:
        - infinity
      resources:
        requests:
          cpu: 4
          memory: 8Gi
'''
            defaultContainer 'astcenc'
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

                    cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../ -DASTCENC_ISA_AVX2=ON ..

                    cov-configure --config ${WORKSPACE}/coverity.conf --template --compiler cc --comptype gcc
                    cov-configure --config ${WORKSPACE}/coverity.conf --template --compiler c++ --comptype g++
                    cov-build --config ${WORKSPACE}/coverity.conf --dir ${WORKSPACE}/intermediate make install
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
            kubernetes {
              yaml '''\
                apiVersion: v1
                kind: Pod
                spec:
                  securityContext:
                    runAsUser: 1000
                    runAsGroup: 1000
                  imagePullSecrets:
                    - name: artifactory-ms-docker
                  containers:
                    - name: astcenc
                      image: mobile-studio--docker.eu-west-1.artifactory.aws.arm.com/astcenc:3.0.0
                      command:
                        - sleep
                      args:
                        - infinity
                      resources:
                        requests:
                          cpu: 4
                          memory: 8Gi
                        limits:
                          cpu: 8
                          memory: 16Gi
              '''.stripIndent()
              defaultContainer 'astcenc'
            }
          }
          stages {
            stage('Clean') {
              steps {
                sh 'git clean -ffdx'
              }
            }
            stage('Build astcenc R x64') {
              steps {
                sh '''
                  export CXX=clang++-9
                  mkdir build_rel
                  cd build_rel
                  cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../ -DASTCENC_ISA_AVX2=ON -DASTCENC_ISA_SSE41=ON -DASTCENC_ISA_SSE2=ON -DASTCENC_PACKAGE=x64 ..
                  make install package -j4
                '''
              }
            }
            stage('Build astcdec R x64') {
              steps {
                sh '''
                  export CXX=clang++-9
                  mkdir build_reldec
                  cd build_reldec
                  cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../ -DASTCENC_ISA_AVX2=ON -DASTCENC_ISA_SSE41=ON -DASTCENC_ISA_SSE2=ON -DASTCENC_DECOMPRESSOR=ON ..
                  make -j4
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
                  python3 ./Test/astc_test_image.py --encoder=all-x86 --test-set Small
                '''
              }
            }
          }
        }
        /* Build for Windows on x86-64 using MSVC ClangCL */
        stage('Windows') {
          agent {
            label 'Windows'
          }
          stages {
            stage('Clean') {
              steps {
                bat 'git clean -ffdx'
              }
            }
            stage('Build R x64') {
              steps {
                bat '''
                  call c:\\progra~2\\micros~1\\2022\\buildtools\\vc\\auxiliary\\build\\vcvars64.bat
                  mkdir build_rel
                  cd build_rel
                  cmake -G "Visual Studio 17 2022" -T ClangCL -DCMAKE_INSTALL_PREFIX=../ -DASTCENC_ISA_AVX2=ON -DASTCENC_ISA_SSE41=ON -DASTCENC_ISA_SSE2=ON -DASTCENC_PACKAGE=x64 ..
                  msbuild astcencoder.sln -property:Configuration=Release
                  msbuild PACKAGE.vcxproj -property:Configuration=Release
                  msbuild INSTALL.vcxproj -property:Configuration=Release
                '''
              }
            }
            stage('Build R Arm64') {
              steps {
                bat '''
                  call c:\\progra~2\\micros~1\\2022\\buildtools\\vc\\auxiliary\\build\\vcvarsall.bat x64_arm64
                  mkdir build_rel_arm64
                  cd build_rel_arm64
                  cmake -G "Visual Studio 17 2022" -A ARM64 -T ClangCL -DCMAKE_INSTALL_PREFIX=../ -DASTCENC_ISA_NEON=ON -DASTCENC_PACKAGE=arm64 ..
                  msbuild astcencoder.sln -property:Configuration=Release
                  msbuild PACKAGE.vcxproj -property:Configuration=Release
                  msbuild INSTALL.vcxproj -property:Configuration=Release
                '''
              }
            }
            stage('Sign') {
              steps {
                dir('sign_tools') {
                    checkout changelog: false,
                             poll: false,
                             scm: [$class: 'GitSCM',
                                   branches: [[name: '*/main']],
                                   doGenerateSubmoduleConfigurations: false,
                                   extensions: [],
                                   submoduleCfg: [],
                                   userRemoteConfigs: [[credentialsId: 'gerrit-jenkins-ssh',
                                                        url: 'ssh://mirror.eu-west-1.gerrit-eu01.aws.arm.com:29418/Hive/shared/signing']]]
                }
                withCredentials([usernamePassword(credentialsId: 'cepe-artifactory-jenkins',
                                                  usernameVariable: 'AF_USER',
                                                  passwordVariable: 'APIKEY')]) {
                    powershell 'C:\\Python311\\python.exe .\\sign_tools\\windows-client-wrapper.py -b $Env:BUILD_NUMBER -t $Env:APIKEY (Get-ChildItem -Filter build_rel\\*.zip)[0].FullName'
                    powershell 'C:\\Python311\\python.exe .\\sign_tools\\windows-client-wrapper.py -b $Env:BUILD_NUMBER -t $Env:APIKEY (Get-ChildItem -Filter build_rel_arm64\\*.zip)[0].FullName'
                }
              }
            }
            stage('Stash') {
              steps {
                dir('build_rel') {
                  stash name: 'astcenc-windows-x64', includes: '*.zip'
                  stash name: 'astcenc-windows-x64-hash', includes: '*.zip.sha256'
                }
                dir('build_rel_arm64') {
                  stash name: 'astcenc-windows-arm64', includes: '*.zip'
                  stash name: 'astcenc-windows-arm64-hash', includes: '*.zip.sha256'
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
            label 'mac && x86_64 && notarizer'
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
                  cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../ -DASTCENC_PACKAGE=universal ..
                  make install package -j4
                '''
              }
            }
            stage('Sign and notarize') {
              environment {
                NOTARIZATION_CREDS = credentials('notarization-account')
              }
              steps {
                dir('build_rel') {
                  sh 'git clone ssh://eu-gerrit-1.euhpc.arm.com:29418/Hive/shared/signing'
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
                  stash name: 'astcenc-macos-universal', includes: '*.zip'
                  stash name: 'astcenc-macos-universal-hash', includes: '*.zip.sha256'
                }
              }
            }
            stage('Test') {
              steps {
                sh '''
                  export PATH=/usr/local/bin:$PATH
                  python3 ./Test/astc_test_image.py --test-set Small --encoder universal
                '''
              }
            }
          }
        }
      }
    }
    stage('Artifactory') {
      agent {
        kubernetes {
          yaml '''
apiVersion: v1
kind: Pod
spec:
  securityContext:
    runAsUser: 1000
    runAsGroup: 1000
  imagePullSecrets:
    - name: artifactory-ms-docker
  containers:
    - name: astcenc
      image: mobile-studio--docker.eu-west-1.artifactory.aws.arm.com/astcenc:3.0.0
      command:
        - sleep
      args:
        - infinity
      resources:
        requests:
          cpu: 1
          memory: 4Gi
'''
          defaultContainer 'astcenc'
        }
      }
      options {
        skipDefaultCheckout true
      }
      stages {
        stage('Unstash') {
          steps {
            dir('upload') {
              unstash 'astcenc-windows-x64-hash'
              unstash 'astcenc-windows-arm64-hash'
              unstash 'astcenc-linux-x64-hash'
              unstash 'astcenc-macos-universal-hash'

              unstash 'astcenc-windows-x64'
              unstash 'astcenc-windows-arm64'
              unstash 'astcenc-linux-x64'
              unstash 'astcenc-macos-universal'

              sh 'cat *.sha256 > release-sha256.txt'
              sh 'rm *.sha256'
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
}
