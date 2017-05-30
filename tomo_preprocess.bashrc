SCRIPTPATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export TOMO_PREPROCESS=$SCRIPTPATH

## PATH

if [ -z ${PATH+x} ]; then
        export PATH=$TOMO_PREPROCESS/bin
else
    	export PATH=$TOMO_PREPROCESS/bin:$PATH
fi

## PYTHONPATH

if [ -z ${PYTHONPATH+x} ]; then
        export PYTHONPATH=$TOMO_PREPROCESS/lib
else
    	export PYTHONPATH=$TOMO_PREPROCESS/lib:$PYTHONPATH
fi
