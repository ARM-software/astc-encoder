def dsgArtifactoryUpload(String sourcePattern, String target) {
    rtBuildInfo (
      // Maximum builds to keep in Artifactory.
      maxBuilds: 10,
      // Also delete the build artifacts when deleting a build.
      deleteBuildArtifacts: true
    )
    rtUpload (
      serverId: 'dsg-artifactory',
      spec: """
        {
          "files": [
            {
              "pattern": "${sourcePattern}",
              "target": "${target}"
            }
          ]
        }
      """
    )
    rtPublishBuildInfo (
      serverId: 'dsg-artifactory'
    )
}

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
        stage('Linux-x86_64') {
          agent {
            docker {
              image 'mobilestudio/astcenc:0.1.0'
              registryUrl 'https://registry.k8s.dsg.arm.com'
              registryCredentialsId 'harbor'
              label 'docker'
              alwaysPull true
            }
          }
          stages {
            stage('Clean workspace') {
              steps {
                sh 'git clean -fdx'
              }
            }
            stage('Linux Build') {
              steps {
                sh '''
                  cd ./Source/
                  make
                '''
              }
            }
            stage('Stash Artefacts') {
              steps {
                dir('Source') {
                  stash name: 'astcenc-linux', includes: 'astcenc'
                }
              }
            }
            stage('Linux Tests') {
              steps {
                sh 'python3 ./Test/astc_test_run.py'
                perfReport(sourceDataFiles:'TestOutput/results.xml')
                junit(testResults: 'TestOutput/results.xml')
              }
            }
          }
        }
        /* Build for Windows on x86_64 */
        stage('Windows-x86_64') {
          agent {
            label 'Windows && x86_64'
          }
          stages {
            stage('Clean workspace') {
              steps {
                bat 'git clean -fdx'
              }
            }
            stage('Windows Release Build') {
              steps {
                bat '''
                  call c:\\progra~2\\micros~1\\2017\\buildtools\\vc\\auxiliary\\build\\vcvars64.bat
                  call msbuild .\\Source\\VS2017\\astcenc.sln /p:Configuration=Release /p:Platform=x64
                '''
              }
            }
            stage('Windows Debug Build') {
              steps {
                bat '''
                  call c:\\progra~2\\micros~1\\2017\\buildtools\\vc\\auxiliary\\build\\vcvars64.bat
                  call msbuild .\\Source\\VS2017\\astcenc.sln /p:Configuration=Debug /p:Platform=x64
                '''
              }
            }
            stage('Stash Artefacts') {
              steps {
                dir('Source\\VS2017\\Release') {
                  stash name: 'astcenc-win-release', includes: 'astcenc.exe'
                }
              }
            }
            stage('Windows Tests') {
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
        stage('MacOS-x86_64') {
          agent {
            label 'mac && x86_64'
          }
          stages {
            stage('Clean workspace') {
              steps {
                sh 'git clean -fdx'
              }
            }
            stage('MacOS Build') {
              steps {
                sh '''
                  cd ./Source/
                  make
                '''
              }
            }
            stage('Stash Artefacts') {
              steps {
                dir('Source') {
                  stash name: 'astcenc-mac', includes: 'astcenc'
                }
              }
            }
            stage('MacOS Tests') {
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
    stage('Upload to Artifactory') {
      agent {
        docker {
          image 'mobilestudio/astcenc:0.1.0'
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
        stage('Unstash artefacts') {
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
        stage('Upload to Artifactory') {
          steps {
            zip zipFile: 'astcenc.zip', dir: 'upload', archive: false
            dsgArtifactoryUpload('*.zip', "astc-encoder/build/${currentBuild.number}/")
          }
        }
      }
    }
  }
}