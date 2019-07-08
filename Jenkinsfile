pipeline {
  agent any
  options {
    buildDiscarder(logRotator(numToKeepStr: '10', artifactNumToKeepStr: '10'))
  }
  tools {
      msbuild 'MSBuild-15.0'
  }
  stages {
    stage('Build') {
      parallel {
        stage('x64 Release') {
          steps {
            script {
              bat(script: "\"${tool 'MSBuild-15.0'}\" .\\Source\\VS2017\\astcenc.sln /p:Configuration=Release /p:Platform=x64", returnStatus: false, returnStdout: false)
            }
          }
        }
        stage('x64 Debug') {
          steps {
            script {
              bat(script: "\"${tool 'MSBuild-15.0'}\" .\\Source\\VS2017\\astcenc.sln /p:Configuration=Debug /p:Platform=x64", returnStatus: false, returnStdout: false)
            }
          }
        }
      }
    }
    stage('Archive') {
      steps {
        archiveArtifacts(artifacts: 'Source\\VS2017\\Release\\astcenc.exe', onlyIfSuccessful: true)
        archiveArtifacts(artifacts: 'Source\\VS2017\\Debug\\astcenc.exe', onlyIfSuccessful: true)
      }
    }
    stage('Test') {
      steps {
        withPythonEnv('Python-3.7.2') {
          bat(script: 'python -m pip install junit_xml')
          bat(script: 'python -m pip install pillow')
          bat(script: 'python Test/runner.py --test-level all --warmup 1 --repeats 5')
        }
      }
    }
  }
  post {
      always {
        recordIssues(enabledForFailure: true, tool: msBuild())
        perfReport(sourceDataFiles:'TestOutput\\results.xml')
        junit(testResults: 'TestOutput\\results.xml')
      }
  }
  triggers {
    pollSCM('H/15 * * * *')
  }
}
