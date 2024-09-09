# ldsc

Build docker file:
docker build --platform linux/amd64 -t ldsc39 .

Run container local:
docker run -it --name ldsc39_container --platform linux/amd64 ldsc39
or
docker run -d -it --platform linux/amd64 --name ldsc39_container --network ldlink_network -p 5000:5000 ldsc39

In the container terminal:(it will activate ldsc enviroment directly)
enter folder: 1kg_eur
cd 1kg_eur

run sample:
python ../ldsc.py --bfile 22 --l2 --ld-wind-cm 1 --out 22

to call endpoint api:
curl http://localhost:5000/ldscore
