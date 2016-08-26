while true; do
  cd
  cd $1
  # NEWMEAS=`find . -type d -print0 | xargs -0 stat -f "%m %N" | sort -rn | head -1 | cut -f2- -d" "`
  NEWMEAS=`find . -type d -print0 | xargs -0 stat -f "%m %N" | sort -rn | head -1 | cut -f2- -d" "`
  echo $NEWMEAS >> $2
  # echo -e $NEWMEAS
  # NEWSAMP=`find $NEWMEAS -name '*.rmg' -type f -print0 | xargs -0 stat -f "%m %N" | sort -rn | head -1 | cut -f2- -d" "`
  cd $NEWMEAS
  find . -name '*.rmg' -type f -print0 | xargs -0 stat -f "%m %N" | sort -rn | head -1 | cut -f2- -d" " >> $2
  # echo $NEWSAMP
  sleep 30
done
