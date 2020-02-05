node('docker') {

    def image
    final String dockerfile = 'jenkins/build.Dockerfile'

    stage('Checkout code') {
        checkout scm
    }

    stage('Build and Push Docker image') {
        final String imageVersion = sh(script: "grep build.environment.version ${dockerfile} | cut -d'\"' -f 2", returnStdout: true).trim()
        println "Building image version ${imageVersion}"
        currentBuild.displayName = imageVersion

        docker.withRegistry('https://registry.k8s.dsg.arm.com', 'harbor') {
            image = docker.build("mobilestudio/astcenc", "- < ${dockerfile}")
            image.push(imageVersion)
            image.push('latest')
        }
    }

    stage('Clean up') {
        sh("docker rmi -f ${image.imageName()}")
    }
}
