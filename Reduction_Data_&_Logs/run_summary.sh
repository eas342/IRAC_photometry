#!/bin/sh

echo
echo
echo 'Would you like to see the list of field names? (y/n)'
read List
echo

if [ "$List" = "y" ]; then
    echo
    echo "The possible field names are:"
    echo "------------------------------"
    echo "1. Date Reduced"
    echo "2. Input Parameters"
    echo "3. Instrument"
    echo "4. File Type"
    echo "5. Target"
    echo "6. Radius Used"
    echo "7. Problem AORs"
    echo "8. Comments"
    echo "Names are case sensitive and must be spelled correctly."
    echo
fi

echo
echo 'Input the field name you would like to see for every run:'
read Name
echo

echo
grep "$Name" *.txt
echo