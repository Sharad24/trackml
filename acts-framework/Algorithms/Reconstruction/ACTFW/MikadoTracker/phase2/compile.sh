sudo docker run -it \
    -v$(pwd):/home/code \
    -v$(pwd)/..:/home/tracker \
    --rm estradevictorantoine/trackml:1.0 \
    /bin/bash -c "cd /home/code; python setup.py build_ext --inplace; exit"
