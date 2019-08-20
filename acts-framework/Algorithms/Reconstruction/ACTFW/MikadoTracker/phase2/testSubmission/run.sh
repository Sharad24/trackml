
INPUT_DATA="/home/sgorbunov/TrackML/data/sample_1000" # CHANGE TO YOUR LOCAL DIRECTORY WITH DATA (hits, cells, and truth csv files)

unzip submission.zip

sudo docker run -i --rm \
	-v $(pwd):/home/code \
	-v $INPUT_DATA:/home/data \
	--cpus=2 --memory=4g \
	estradevictorantoine/trackml:1.0 \
	/bin/sh -c "cd /home/code; python main.py $*"
