pipeline {
  agent any
  triggers {
    pollSCM('H/15 * * * *')
  }
  stages {
    stage('Build Release') {
      steps {
        echo 'Build test'
        script {
          def MSBuild = tool(name: 'MSBuild', type: 'hudson.plugins.msbuild.MsBuildInstallation')
          bat(script: "${MSBuild} .\\Source\\win32-2017\\astcenc\\astcenc.sln /p:Configuration=Release /p:Platform=x64", returnStatus: true, returnStdout: true)
        }
      }
    }
  }
}
