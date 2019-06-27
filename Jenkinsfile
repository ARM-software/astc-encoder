pipeline {
  agent any
  stages {
    stage('Build Release') {
      steps {
        echo 'Build test'
        bat(script: '\\"${tool \'MSBuild\'}\\" .\\Source\\win32-2017\\astcenc\\astcenc.sln /p:Configuration=Release /p:Platform=x64', returnStatus: true, returnStdout: true)
      }
    }
  }
}