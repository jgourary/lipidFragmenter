package main

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"os"
	filepath2 "path/filepath"
	"strconv"
	"strings"
)

func generateLibrary(inPath string, outPath string, outDir string, limit int) {

	_ = os.MkdirAll(outDir, 0755)

	// open file
	file, err := os.Open(inPath)
	if err != nil {
		fmt.Println("Failed to open molecule file: " + inPath)
		log.Fatal(err)
	}

	// open file
	outFile, err := os.Create(outPath)
	if err != nil {
		fmt.Println("Failed to open molecule file: " + outPath)
		log.Fatal(err)
	}

	// Initialize scanner
	scanner := bufio.NewScanner(file)
	// iterate over all other lines
	for scanner.Scan() {
		// get next line
		line := scanner.Text()
		// split by whitespace
		tokens := strings.Fields(line)
		if len(tokens) > 2 {
			tokens := strings.Fields(line)

			smiles := tokens[0]
			val, _ := strconv.Atoi(tokens[1])
			if val < limit {
				break
			}
			path := tokens[2]
			name := filepath2.Base(path)
			baseName := strings.Split(name, ".")[0]
			dir := filepath2.Join(outDir, baseName)
			_ = os.Mkdir(dir, 0755)
			out := filepath2.Join(dir, name)
			_, err = copyFile(path, out)

			_, _ = outFile.WriteString(smiles + "\t" + out + "\n")
		}

	}

}

func copyFile(src, dst string) (int64, error) {
	sourceFileStat, err := os.Stat(src)
	if err != nil {
		return 0, err
	}

	if !sourceFileStat.Mode().IsRegular() {
		return 0, fmt.Errorf("%s is not a regular file", src)
	}

	source, err := os.Open(src)
	if err != nil {
		return 0, err
	}
	defer source.Close()

	destination, err := os.Create(dst)
	if err != nil {
		return 0, err
	}
	defer destination.Close()
	nBytes, err := io.Copy(destination, source)
	return nBytes, err
}
