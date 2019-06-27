pipeline {
  agent any
  stages {
    stage('Build') {
      parallel {
        stage('x64 Release') {
          steps {
            echo 'Build test'
            script {
              def MSBuild = tool(name: 'MSBuild', type: 'hudson.plugins.msbuild.MsBuildInstallation')
              bat(script: "\"${MSBuild}\" .\\Source\\win32-2017\\astcenc\\astcenc.sln /p:Configuration=Release /p:Platform=x64", returnStatus: true, returnStdout: true)
            }

          }
        }
        stage('x64 Debug') {
          steps {
            script {
              def MSBuild = tool(name: 'MSBuild', type: 'hudson.plugins.msbuild.MsBuildInstallation')
              bat(script: "\"${MSBuild}\" .\\Source\\win32-2017\\astcenc\\astcenc.sln /p:Configuration=Debug /p:Platform=x64", returnStatus: true, returnStdout: true)
            }

          }
        }
        stage('x86 Release') {
          steps {
            script {
              def MSBuild = tool(name: 'MSBuild', type: 'hudson.plugins.msbuild.MsBuildInstallation')
              bat(script: "\"${MSBuild}\" .\\Source\\win32-2017\\astcenc\\astcenc.sln /p:Configuration=Release /p:Platform=Win32", returnStatus: true, returnStdout: true)
            }

          }
        }
        stage('x86 Debug') {
          steps {
            script {
              def MSBuild = tool(name: 'MSBuild', type: 'hudson.plugins.msbuild.MsBuildInstallation')
              bat(script: "\"${MSBuild}\" .\\Source\\win32-2017\\astcenc\\astcenc.sln /p:Configuration=Debug /p:Platform=Win32", returnStatus: true, returnStdout: true)
            }

          }
        }
      }
    }
    stage('Archive') {
      steps {
        archiveArtifacts(onlyIfSuccessful: true, artifacts: 'Source\\win32-2017\\astcenc\\Win32\\Release\\astcenc.exe')
        archiveArtifacts(artifacts: 'Source\\win32-2017\\astcenc\\Win32\\Debug\\astcenc.exe', onlyIfSuccessful: true)
        archiveArtifacts(artifacts: 'Source\\win32-2017\\astcenc\\x64\\Release\\astcenc.exe', onlyIfSuccessful: true)
        archiveArtifacts 'Source\\win32-2017\\astcenc\\x64\\Debug\\astcenc.exe'
      }
    }
  }
  triggers {
    pollSCM('H/15 * * * *')
  }
}