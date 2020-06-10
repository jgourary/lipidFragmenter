package main

import (
	"bufio"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	filepath2 "path/filepath"
	"strconv"
	"strings"
)

func smiProcessor(dir string, fragmentsDir string) {

	fileInfo, err := ioutil.ReadDir(fragmentsDir)
	if err != nil {
		fmt.Println("failed to read directory: " + fragmentsDir)
		log.Fatal(err)
	}
	smiDict := make(map[string]int)

	for i := 0; i < len(fileInfo); i++ {
		// if item is a Dir (as it should be unless the end user tampered with the directory manually...)
		if fileInfo[i].IsDir() {
			fragSubDir := filepath2.Join(fragmentsDir, fileInfo[i].Name())

			fragSubDirInfo, err := ioutil.ReadDir(fragSubDir)
			if err != nil {
				fmt.Println("failed to read directory: " + fragSubDir)
				log.Fatal(err)
			}
			for j := 0; j < len(fragSubDirInfo); j++ {
				fileName := fragSubDirInfo[j].Name()
				if filepath2.Ext(fileName) == ".can" {

					// open file
					smiFilePath := filepath2.Join(fragSubDir,fileName)
					file, err := os.Open(smiFilePath)
					if err != nil {
						fmt.Println("Failed to open molecule file: " + smiFilePath)
						log.Fatal(err)
					}
					// Initialize scanner
					scanner := bufio.NewScanner(file)
					// ignore first line
					scanner.Scan()
					line := scanner.Text()
					fields := strings.Fields(line)
					if len(fields) > 0 {

						smiString := fields[0]
						isAlkane := true
						alkaneRunes := "CcHh()[]1234567890@/="
						backslashRune := []rune("\\")[0]
						alkaneRunes2 := []rune(alkaneRunes)
						alkaneRunes2 = append(alkaneRunes2, backslashRune)

						for _, rune := range smiString {
							runeIsAlkaneValid := false
							for _, alkaneRune := range alkaneRunes2 {
								if rune == alkaneRune {
									runeIsAlkaneValid = true
									break
								}
							}
							if !runeIsAlkaneValid {
								isAlkane = false
								break
							}
						}

						if isAlkane == false {
							// if group not in map already
							if _, ok := smiDict[fields[0]]; !ok {
								// add group to map
								smiDict[fields[0]] = 1
							} else {
								// append atom to group in map
								smiDict[fields[0]] += 1
							}
						}
					}
				}
			}
		}
	}

	thisPath := filepath2.Join(dir, "fragments.txt")
	// fmt.Println(thisPath)
	thisFile, err := os.Create(thisPath)
	if err != nil {
		fmt.Println("Failed to create new file: " + thisPath)
		log.Fatal(err)
	}
	
	keys := make([]string, len(smiDict))
	vals := make([]int, len(smiDict))

	i := 0
	for key, val := range smiDict {
		keys[i] = key
		vals[i] = val
		i++
	}
	vals, keys = qsort(vals, keys)
	keys = reverseString(keys)
	vals = reverseInt(vals)

	for i := 0; i < len(keys); i++ {
		_, err = thisFile.WriteString(keys[i] + "\t" + strconv.Itoa(vals[i]) + "\n")
		if err != nil {
			fmt.Println("Failed to write line to file: " + thisPath)
			log.Fatal(err)
		}
	}
	fmt.Println("Number of Unique Fragments = " + strconv.Itoa(len(keys)) )
	numFrags := 0
	for _, val := range vals { numFrags += val }
	fmt.Println("Total Number of Fragments in " +  strconv.Itoa(len(fileInfo)) + " molecules = " + strconv.Itoa(numFrags))

	numFragsArray := []int{10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200}
	for i := 0 ; i < len(numFragsArray); i++ {
		FragsCovered := 0
		for j := 0 ; j < numFragsArray[i]; j++ {
			FragsCovered += vals[j]
		}
		covPerc := 100.0 * float64(FragsCovered) / float64(numFrags)
		covString := fmt.Sprintf("%.2f",covPerc)
		fmt.Println("Top " + strconv.Itoa(numFragsArray[i]) + " fragments represent " + covString + "% of total fragment occurrences")
	}

	//smiCoverage(fragmentsDir,keys)
}


func smiCoverage(fragmentsDir string, keys []string) {



	numFrags := []int{10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200}

	thisMap := make(map[string]int)

	for i := 0 ; i < len(numFrags); i++ {
		var min int
		if i == 0 {
			min = 0
		} else {
			min = numFrags[i-1]
		}
		for j := min ; j < numFrags[i]; j++ {
			thisMap[keys[j]] = 0
		}
		fmt.Println(len(thisMap))
		coverage := getCoverage(fragmentsDir,thisMap)
		fmt.Println("Top " + strconv.Itoa(numFrags[i]) + " fragments: " + fmt.Sprintf("%.2f",100*coverage) + "%")
	}

}

func getCoverage(fragmentsDir string, thisMap map[string]int) float64 {
	fileInfo, err := ioutil.ReadDir(fragmentsDir)
	if err != nil {
		fmt.Println("failed to read directory: " + fragmentsDir)
		log.Fatal(err)
	}

	numNotInSet := 0
	totalNum := len(fileInfo)
	for i := 0; i < len(fileInfo); i++ {
		// if item is a Dir (as it should be unless the end user tampered with the directory manually...)
		if fileInfo[i].IsDir() {
			fragSubDir := filepath2.Join(fragmentsDir, fileInfo[i].Name())

			fragSubDirInfo, err := ioutil.ReadDir(fragSubDir)
			if err != nil {
				fmt.Println("failed to read directory: " + fragSubDir)
				log.Fatal(err)
			}
			for j := 0; j < len(fragSubDirInfo); j++ {
				fileName := fragSubDirInfo[j].Name()
				if filepath2.Ext(fileName) == ".smi" {

					// open file
					smiFilePath := filepath2.Join(fragSubDir, fileName)
					file, err := os.Open(smiFilePath)
					if err != nil {
						fmt.Println("Failed to open molecule file: " + smiFilePath)
						log.Fatal(err)
					}
					// Initialize scanner
					scanner := bufio.NewScanner(file)
					// ignore first line
					scanner.Scan()
					line := scanner.Text()
					fields := strings.Fields(line)
					if len(fields) > 0 {
						smiString := fields[0]
						// if group not in map already
						if _, ok := thisMap[smiString]; !ok {
							// add group to map
							numNotInSet++
							break
						}

					}
				}
			}
		}
	}
	return 1.0 - (float64(numNotInSet) / float64(totalNum))

}





