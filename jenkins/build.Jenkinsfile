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
        /* Build for Linux on x86_64 */
        stage('Linux') {
          agent {
            docker {
              image 'mobilestudio/astcenc:0.2.0'
              registryUrl 'https://registry.k8s.dsg.arm.com'
              registryCredentialsId 'harbor'
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
            stage('Build') {
              steps {
                sh '''
                  cd ./Source/
                  make CXX=clang++ VEC=avx2
                '''
              }
            }
            stage('Stash') {
              steps {
                dir('Source') {
                  stash name: 'astcenc-linux', includes: 'astcenc-*'
                }
              }
            }
            stage('Test') {
              steps {
                sh 'python3 ./Test/astc_test_run.py'
                perfReport(sourceDataFiles:'TestOutput/results.xml')
                junit(testResults: 'TestOutput/results.xml')
              }
            }
          }
        }
        /* Build for Windows on x86_64 */
        stage('Windows') {
          agent {
            label 'Windows && x86_64'
          }
          stages {
            stage('Clean') {
              steps {
                bat 'git clean -fdx'
              }
            }
            stage('Build R') {
              steps {
                bat '''
                  call c:\\progra~2\\micros~1\\2019\\buildtools\\vc\\auxiliary\\build\\vcvars64.bat
                  call msbuild .\\Source\\VS2019\\astcenc-avx2.vcxproj /p:Configuration=Release /p:Platform=x64
                '''
              }
            }
            stage('Build D') {
              steps {
                bat '''
                  call c:\\progra~2\\micros~1\\2019\\buildtools\\vc\\auxiliary\\build\\vcvars64.bat
                  call msbuild .\\Source\\VS2019\\astcenc-avx2.vcxproj /p:Configuration=Debug /p:Platform=x64
                '''
              }
            }
            stage('Stash') {
              steps {
                dir('Source\\VS2019\\astcenc-avx2-Release') {
                  stash name: 'astcenc-win-release', includes: 'astcenc-avx2.exe'
                }
              }
            }
            stage('Test') {
              steps {
                bat '''
                  set Path=c:\\Python38;c:\\Python38\\Scripts;%Path%
                  call python ./Test/astc_test_run.py
                '''
                perfReport(sourceDataFiles:'TestOutput\\results.xml')
                junit(testResults: 'TestOutput\\results.xml')
              }
            }
          }
        }
        /* Build for Mac on x86_64 */
        stage('macOS') {
          agent {
            label 'mac && x86_64'
          }
          stages {
            stage('Clean') {
              steps {
                sh 'git clean -fdx'
              }
            }
            stage('Build') {
              steps {
                sh '''
                  cd ./Source/
                  make VEC=avx2
                '''
              }
            }
            stage('Stash') {
              steps {
                dir('Source') {
                  stash name: 'astcenc-mac', includes: 'astcenc-*'
                }
              }
            }
            stage('Test') {
              steps {
                sh '''
                  export PATH=$PATH:/usr/local/bin
                  python3 ./Test/astc_test_run.py
                '''
                perfReport(sourceDataFiles:'TestOutput/results.xml')
                junit(testResults: 'TestOutput/results.xml')
              }
            }
          }
        }
      }
    }
    stage('Artifactory') {
      agent {
        docker {
          image 'mobilestudio/astcenc:0.2.0'
          registryUrl 'https://registry.k8s.dsg.arm.com'
          registryCredentialsId 'harbor'
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
            deleteDir()
            dir('upload/Linux-x86_64') {
              unstash 'astcenc-linux'
            }
            dir('upload/Windows-x86_64') {
              unstash 'astcenc-win-release'
            }
            dir('upload/MacOS-x86_64') {
              unstash 'astcenc-mac'
            }
          }
        }
        stage('Upload') {
          steps {
            zip zipFile: 'astcenc.zip', dir: 'upload', archive: false
            dsgArtifactoryUpload('*.zip')
            dsgArtifactoryPromote()
          }
        }
      }
    }
  }
}