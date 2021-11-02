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

func fragSelector(dir string, singleFragmentsDir string, doubleFragmentsDir string, uniqueSFDir string, uniqueDFDir string) {

	singleFragRankedKeys, singleFragRankedVals, singleFragStringToFragLocations, isSFhydrocarbon := countFrags(singleFragmentsDir)
	outPath := filepath2.Join(dir,"top_single_fragments.txt")
	outPathHC := filepath2.Join(dir,"top_single_fragments_HC.txt")
	writeTopFrags(outPath, outPathHC, singleFragRankedKeys, singleFragRankedVals, singleFragStringToFragLocations, isSFhydrocarbon)

	doubleFragRankedKeys, doubleFragRankedVals, doubleFragStringToFragLocations, isDFhydrocarbon := countFrags(doubleFragmentsDir)
	outPath = filepath2.Join(dir,"top_double_fragments.txt")
	outPathHC = filepath2.Join(dir,"top_double_fragments_HC.txt")
	writeTopFrags(outPath, outPathHC, doubleFragRankedKeys, doubleFragRankedVals, doubleFragStringToFragLocations, isDFhydrocarbon)

	_ = os.MkdirAll(uniqueSFDir, 0755)
	_ = os.MkdirAll(uniqueDFDir, 0755)
	for i := 0; i < len(singleFragRankedKeys); i++ {
		thisPath := filepath2.Join(uniqueSFDir, "unique_single_" + strconv.Itoa(i) + ".info")
		// fmt.Println(thisPath)
		thisFile, err := os.Create(thisPath)
		if err != nil {
			fmt.Println("Failed to create new unique single fragment file: " + thisPath)
			log.Fatal(err)
		}
		var HCstatus string
		if isSFhydrocarbon[singleFragRankedKeys[i]] == true { HCstatus = "hydrocarbon"} else { HCstatus = "non-hydrocarbon"}
		_, err = thisFile.WriteString(singleFragRankedKeys[i] + "\t" + HCstatus + "\t" + strconv.Itoa(singleFragRankedVals[i]) + "\n")
		for j:=0; j < len(singleFragStringToFragLocations[singleFragRankedKeys[i]]); j++ {
			_, err = thisFile.WriteString(singleFragStringToFragLocations[singleFragRankedKeys[i]][j] + "\n")
		}
	}

	for i := 0; i < len(doubleFragRankedKeys); i++ {
		thisPath := filepath2.Join(uniqueDFDir, "unique_double_" + strconv.Itoa(i) + ".info")
		// fmt.Println(thisPath)
		thisFile, err := os.Create(thisPath)
		if err != nil {
			fmt.Println("Failed to create new double fragment file: " + thisPath)
			log.Fatal(err)
		}
		var HCstatus string
		if isDFhydrocarbon[doubleFragRankedKeys[i]] == true { HCstatus = "hydrocarbon"} else { HCstatus = "non-hydrocarbon"}
		_, err = thisFile.WriteString(doubleFragRankedKeys[i] + "\t" + HCstatus + "\t" + strconv.Itoa(doubleFragRankedVals[i]) + "\n")
		for j:=0; j < len(doubleFragStringToFragLocations[doubleFragRankedKeys[i]]); j++ {
			_, err = thisFile.WriteString(doubleFragStringToFragLocations[doubleFragRankedKeys[i]][j] + "\n")
		}
	}
}
func countFrags(fragmentsDir string) ([]string, []int, map[string][]string, map[string]bool) {

	fragDirFileInfo, err := ioutil.ReadDir(fragmentsDir)
	if err != nil {
		fmt.Println("failed to read directory: " + fragmentsDir)
		log.Fatal(err)
	}

	fragStringToFragCount := make(map[string]int)
	fragStringToFragLocations := make(map[string][]string)
	isFragHydrocarbon := make(map[string]bool)



	fmt.Println("Identifying all non-alkane fragments...")

	for i := 0; i < len(fragDirFileInfo); i++ {
		fileName := fragDirFileInfo[i].Name()
		if filepath2.Ext(fileName) == ".can" {

			canFilePath := filepath2.Join(fragmentsDir, fileName)
			txyzFilePath := strings.Split(canFilePath,".")[0] + ".txyz"


			smiString := getSMIString(canFilePath)
			smiString = makeSMILESUnique(smiString)
			isHydrocarbon := isSMILESHydrocarbon(smiString)
			
			// if smiles string not in map already
			if _, ok := fragStringToFragCount[smiString]; !ok {
				// add string to maps
				fragStringToFragCount[smiString] = 1
				fragStringToFragLocations[smiString] = []string{txyzFilePath}
				isFragHydrocarbon[smiString] = isHydrocarbon
			} else {
				// amend entry of string in maps
				fragStringToFragCount[smiString] += 1
				fragStringToFragLocations[smiString] = append(fragStringToFragLocations[smiString], txyzFilePath)
			}
			
		}
			
		
	}

	fmt.Println("Sorting fragments...")

	keys := make([]string, len(fragStringToFragCount))
	vals := make([]int, len(fragStringToFragCount))

	i := 0
	for key, val := range fragStringToFragCount {
		keys[i] = key
		vals[i] = val
		i++
	}
	vals, keys = qsort(vals, keys)
	keys = reverseString(keys)
	vals = reverseInt(vals)

	return keys, vals, fragStringToFragLocations, isFragHydrocarbon
}

func isSMILESHydrocarbon(smiString string) bool {

	hydrocarbonRunes := "CcHh()[]1234567890@/="
	backslashRune := []rune("\\")[0]
	hydrocarbonRunes2 := []rune(hydrocarbonRunes)
	hydrocarbonRunes2 = append(hydrocarbonRunes2, backslashRune)

	isHydrocarbon := true

	for _, thisRune := range smiString {
		runeIsHydrocarbonValid := false
		for _, hydrocarbonRune := range hydrocarbonRunes2 {
			if thisRune == hydrocarbonRune {
				runeIsHydrocarbonValid = true
				break
			}
		}
		if !runeIsHydrocarbonValid {
			isHydrocarbon = false
			break
		}
	}
	return isHydrocarbon
}

func makeSMILESUnique(smiString string) string {
	searchPhrase := "([O])O"
	replacePhrase := "(O)[O]"

	return strings.ReplaceAll(smiString, searchPhrase, replacePhrase)
}

// retrieves the canonical smiles string from a .can file
func getSMIString(smiFilePath string) string {
	// open file
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
		return smiString
	} else {
		return "null"
	}
}

func writeTopFrags(outPath string, outPathHC string, topKeys []string, topVals []int, keyToLocns map[string][]string, isHC map[string]bool) {

	outFile, _ := os.Create(outPath)
	outFileHC, _ := os.Create(outPathHC)

	for i := 0; i < len(topKeys); i++ {

		if isHC[topKeys[i]] == true {
			_, _ = outFileHC.WriteString(topKeys[i] + "\t" + strconv.Itoa(topVals[i]) + "\t" + keyToLocns[topKeys[i]][0] + "\n")
		} else {
			_, _ = outFile.WriteString(topKeys[i] + "\t" + strconv.Itoa(topVals[i]) + "\t" + keyToLocns[topKeys[i]][0] + "\n")
		}

	}
}