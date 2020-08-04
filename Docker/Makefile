ceo: Dockerfile
	cd ../; git archive --format=tar.gz HEAD > Docker/ceo.tar.gz
	sudo docker build -t gmto.im/ceo:devel_cuda10.1 .

cluster:
	eksctl create cluster --name DOS --region us-west-2 --nodegroup-name DOSlet --node-type p3.2xlarge --nodes 1 --nodes-min 0 --nodes-max 4 --node-ami ami-0ca5998dc2c88e64b
	kubectl apply -f https://raw.githubusercontent.com/NVIDIA/k8s-device-plugin/1.0.0-beta/nvidia-device-plugin.yml

test:
	kubectl apply -f ceo-test.yaml

clean:
	eksctl delete cluster --name DOS --region us-west-2
