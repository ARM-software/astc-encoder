pipeline {
  agent any
  stages {
    stage('Build') {
      parallel {
        stage('Build Release') {
          steps {
            echo 'Build test'
            script {
              def MSBuild = tool(name: 'MSBuild', type: 'hudson.plugins.msbuild.MsBuildInstallation')
              bat(script: "\"${MSBuild}\" .\\Source\\win32-2017\\astcenc\\astcenc.sln /p:Configuration=Release /p:Platform=x64", returnStatus: true, returnStdout: true)
            }

          }
        }
        stage('Build Debug') {
          steps {
            script {
              def MSBuild = tool(name: 'MSBuild', type: 'hudson.plugins.msbuild.MsBuildInstallation')
              bat(script: "\"${MSBuild}\" .\\Source\\win32-2017\\astcenc\\astcenc.sln /p:Configuration=Debug /p:Platform=x64", returnStatus: true, returnStdout: true)
            }

          }
        }
      }
    }
    stage('Archive') {
      steps {
        archiveArtifacts(onlyIfSuccessful: true, artifacts: 'Source\\win32-2017\\astcenc\\x64\\Release\\astcenc.exe,Source\\win32-2017\\astcenc\\x64\\Debug\\astcenc.exe,')
      }
    }
  }
  triggers {
    pollSCM('H/15 * * * *')
  }
}