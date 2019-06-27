pipeline {
  agent any
  stages {
    stage('Build Release') {
      steps {
        echo 'Build test'
        sh '''msbuild .\\Source\\win32-2017\\astcenc\\astcenc.sln \\
    /p:Configuration=Release \\
    /p:Platform=x64'''
      }
    }
  }
}