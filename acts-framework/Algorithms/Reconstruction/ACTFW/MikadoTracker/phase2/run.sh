
INPUT_DATA="/home/sgorbunov/TrackML/data/sample_1000" # CHANGE TO YOUR LOCAL DIRECTORY WITH DATA (hits, cells, and truth csv files)
cp ../geoLayerSizes.txt .
cp ../geoLayerModules.txt .
cp ../geoLayerField.txt .
cp ../cuts.txt .


sudo docker run -i --rm \
	-v $(pwd):/home/code \
	-v $INPUT_DATA:/home/data \
	--cpus=2 --memory=4g \
	estradevictorantoine/trackml:1.0 \
	/bin/sh -c "cd /home/code; python main.py $*"
#	/bin/sh -c "cd /home/code; export LD_LIBRARY_PATH=/home/code/:$LD_LIBRARY_PATH; python main.py $*"
