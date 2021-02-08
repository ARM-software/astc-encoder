pipeline {
  agent {
    kubernetes {
      yaml '''
apiVersion: v1
kind: Pod
spec:
  securityContext:
    privileged: true
  imagePullSecrets:
    - name: artifactory-ms-docker
  containers:
    - name: dind
      image: mobile-studio--docker.eu-west-1.artifactory.aws.arm.com/docker:dind
      command:
        - sleep
      args:
        - infinity
      resources:
        requests:
          cpu: 2
          memory: 4Gi
'''
            defaultContainer 'dind'
    }
  }
  environment {
    ARTIFACTORY_CREDENTIALS = credentials('cepe-artifactory-jenkins')
  }
  options {
    ansiColor('xterm')
    timestamps()
  }
  stages {
    stage('Build and Push Image') {
      steps {
        sh '''
          apk add --no-cache bash
          chmod u+x ./jenkins/build-image.sh
          ./jenkins/build-image.sh push
        '''
      }
    }
  }
}
