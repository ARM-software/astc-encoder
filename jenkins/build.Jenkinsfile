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
            stage('Build') {
              steps {
                sh '''
                  cd ./Source/
                  make
                '''
              }
            }
            stage('Archive Artefacts') {
              steps {
                archiveArtifacts(artifacts: './Source/astcenc', onlyIfSuccessful: true)
              }
            }
            stage('Test') {
              steps {
                sh 'python3 ./Test/astc_test_run.py'
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
            stage('Release') {
              steps {
                bat '''
                  call c:\\progra~2\\micros~1\\2017\\buildtools\\vc\\auxiliary\\build\\vcvars64.bat
                  call msbuild .\\Source\\VS2017\\astcenc.sln /p:Configuration=Release /p:Platform=x64
                '''
              }
            }
            stage('Debug') {
              steps {
                bat '''
                  call c:\\progra~2\\micros~1\\2017\\buildtools\\vc\\auxiliary\\build\\vcvars64.bat
                  call msbuild .\\Source\\VS2017\\astcenc.sln /p:Configuration=Debug /p:Platform=x64
                '''
              }
            }
            stage('Archive Artefacts') {
              steps {
                archiveArtifacts(artifacts: 'Source\\VS2017\\Release\\astcenc.exe', onlyIfSuccessful: true)
                archiveArtifacts(artifacts: 'Source\\VS2017\\Debug\\astcenc.exe', onlyIfSuccessful: true)
              }
            }
            stage('Test') {
              steps {
                bat '''
                  set Path=c:\\Python38;c:\\Python38\\Scripts;%Path%
                  call python ./Test/astc_test_run.py
                '''
              }
            }
          }
        }
        stage('MacOS-x86_64') {
          agent {
            label 'mac && x86_64'
          }
          stages {
            stage('Build') {
              steps {
                sh '''
                  cd ./Source/
                  make
                '''
              }
            }
            stage('Archive Artefacts') {
              steps {
                archiveArtifacts(artifacts: './Source/astcenc', onlyIfSuccessful: true)
              }
            }
            stage('Test') {
              steps {
                sh 'python3 ./Test/astc_test_run.py'
              }
            }
          }
        }
      }
    }
  }
}