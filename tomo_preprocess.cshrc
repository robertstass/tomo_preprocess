set SCRIPT = `/usr/sbin/lsof +p $$ | grep -oE /.\*tomo_preprocess.cshrc`
set SCRIPTPATH = `dirname $SCRIPT`
setenv TOMO_PREPROCESS $SCRIPTPATH


### Paths ###

if ( $?PATH ) then
        setenv PATH $TOMO_PREPROCESS/bin:$PATH
else
        setenv PATH $TOMO_PREPROCESS/bin
endif

### Python path
if ( $?PYTHONPATH ) then
	setenv PYTHONPATH $TOMO_PREPROCESS/python_dependencies:$PYTHONPATH
else
	setenv PYTHONPATH $TOMO_PREPROCESS/python_dependencies
endif





