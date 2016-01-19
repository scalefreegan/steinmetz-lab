30 03 01 */3 * /path/to/script

find ./ -user brooks -print0 | xargs -0 du -hc | tail -n1

find ./ -type f -size +20M -user brooks 2> /dev/null  -print0 | xargs -0 du -hc | tail -n1

find ./ -type f -size +20M -user brooks -execdir ls -lh {} \; 2> /dev/null | awk '{ print $NF ": " $5 }' | sort -nk 2,2

echo "hi" | mailx -v -s "Space usage on /g/steinmetz" thetruemissinglink@gmail.com
