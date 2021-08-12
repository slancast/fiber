#! bash
#rename files from list
#found on: https://unix.stackexchange.com/questions/152656/rename-a-list-of-files-according-to-a-text-file

for file in *; do mv "${file}" "${file//_2/}"; done

for file in *.txt; do read -r line;  mv -v "${file}" "${line}";  done < ../filenames_lines.csv

#list is text file with list of names
#for some reason this command adds a question mark. Not sure why, but to remove it:

for filename in *
do 
    if [ "$filename" == *"?"* ] 
    then
        mv "$filename" "$(echo $filename | tr '?' '-')" 
    fi
done

