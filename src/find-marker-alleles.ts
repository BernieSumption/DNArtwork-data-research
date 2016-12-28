/// <reference path="../typings/index.d.ts" />

/**
 * Parse a collection of gzip'd hapmap allele frequency files and output a
 * list of alleles that serve as relatively rare markers, i.e. they are within
 * a defined range of frequencies in every popultion.
 */

const MIN_MARKER_FREQ = 0.02;
const MAX_MARKER_FREQ = 0.10;
const MARKERS_PER_CHROMASOME = 200;
const POPULATIONS = ["ASW", "CEU", "CHB", "CHD", "GIH", "JPT", "LWK", "MEX", "MKK", "TSI", "YRI"];
const CHROMASOMES = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"];
// const POPULATIONS = ["ASW", "CEU"];
// const CHROMASOMES = ["3"];
const FREQS_FILE_PATTERN = "./__tmp_hapmap-freqs/allele_freqs_chr{CHR}_{POP}_r28_nr.b36_fwd.txt.gz";
const TESTED_MARKERS_FILE = "./marker-list/markers-list-intersection.txt";

import _ = require("underscore");
import fs = require("fs");
import zlib = require("zlib");
import readline = require("readline");

let analyser: Promise<any>[] = [];

let CHROMASOME: string;
let MARKERS: { [slug: string]: Marker };
let SCRIPT_RESULTS: ChromosomeMarkerList[] = [];

let chrIndex = 0;
let testedMarkersIds: { [rsId: string]: boolean };

recursiveAnalyseChromasomes()
  .then(() => logNotice("DONE!"));

function recursiveAnalyseChromasomes(): Promise<{}> {
  CHROMASOME = CHROMASOMES[chrIndex++];
  if (CHROMASOME) {
    logNotice(`Loading tested markers`);
    loadTestedMarkers();
    logNotice(`Starting analysis of chromasome ${CHROMASOME}`);
    MARKERS = {};
    return chromasomeAnalyser(CHROMASOME)
      .then(collectChromasomeResults)
      .then(recursiveAnalyseChromasomes);
  } else {
    printResults();
  }
}

function loadTestedMarkers() {
  testedMarkersIds = {};
  for (let rsid of fs.readFileSync(TESTED_MARKERS_FILE, "utf8").split("\n")) {
    testedMarkersIds[rsid.trim()] = true;
  }
}

function chromasomeAnalyser(chromasome: string) {
  return Promise.all(POPULATIONS.map(population => populationChromasomeAnalyser(population, chromasome)));
}

function populationChromasomeAnalyser(population: string, chromasome: string) {
  return new Promise((resolve, reject) => {
    let filePath = FREQS_FILE_PATTERN.replace("{POP}", population).replace("{CHR}", chromasome);
    let firstLine = true;
    readline.createInterface({
      input: fs.createReadStream(filePath).pipe(zlib.createGunzip())
    })
      .on("line", (line: string) => {
        if (firstLine) {
          if (line !== "rs# chrom pos strand build center protLSID assayLSID panelLSID QC_code refallele refallele_freq refallele_count otherallele otherallele_freq otherallele_count totalcount") {
            fatalError("Bad header format");
          }
          firstLine = false;
        } else {
          let cells = line.split(" ");
          // logNotice("rsid = " + cells[HapmapFields.RSID] + " tested = " + testedMarkersIds[cells[HapmapFields.RSID]]);
          if (testedMarkersIds[cells[HapmapFields.RSID]]) {
            processMarker(cells[HapmapFields.RSID], cells[HapmapFields.REFALLELE], parseFloat(cells[HapmapFields.REFALLELE_FREQ]), chromasome);
            processMarker(cells[HapmapFields.RSID], cells[HapmapFields.OTHERALLELE], parseFloat(cells[HapmapFields.OTHERALLELE_FREQ]), chromasome);
          }
        }
      })
      .on("error", reject)
      .on("close", () => {
        logNotice(`   Finished processing ${chromasome}>${population}`);
        resolve();
      });
  });

  function processMarker(rsid: string, allele: string, freq: number, chr: string) {
    let slug = `${rsid}/${allele}`;
    let marker = MARKERS[slug];
    if (!marker) {
      marker = MARKERS[slug] = {
        slug: slug,
        chr: chr,
        rsid: rsid,
        allele: allele,
        populations: 0,
        populationFreqs: [],
        minFreq: 1,
        avgFreq: 0,
        maxFreq: 0
      };
    }
    marker.populationFreqs[marker.populations] = freq;
    ++marker.populations;
    let prop = 1 / marker.populations;
    marker.avgFreq *= 1 - prop;
    marker.avgFreq += freq * prop;
    if (marker.minFreq > freq) {
      marker.minFreq = freq;
    }
    if (marker.maxFreq < freq) {
      marker.maxFreq = freq;
    }
  }
}

function collectChromasomeResults() {
  let populationCount = POPULATIONS.length;
  let acceptableMarkers: Marker[] = [];
  for (let slug in MARKERS) {
    let marker = MARKERS[slug];
    if (marker.minFreq > MIN_MARKER_FREQ /*&& marker.maxFreq < MAX_MARKER_FREQ */ && marker.populations === populationCount) {
      acceptableMarkers.push(marker);
    }
  }
  if (acceptableMarkers.length < MARKERS_PER_CHROMASOME) {
    throw new Error(`Not enough acceptable markers on chr ${CHROMASOME} (got ${acceptableMarkers.length}, required ${MARKERS_PER_CHROMASOME})`);
  }
  acceptableMarkers = _.sortBy(acceptableMarkers, m => m.maxFreq);
  let originalCount = acceptableMarkers.length;

  let markersToDisplay = acceptableMarkers;
  if (process.argv.indexOf("--all") === -1) {
    // markersToDisplay = acceptableMarkers.slice(0, MARKERS_PER_CHROMASOME/2).concat(acceptableMarkers.slice(acceptableMarkers.length - MARKERS_PER_CHROMASOME/2));
    markersToDisplay = acceptableMarkers.slice(0, MARKERS_PER_CHROMASOME);
  }

  if (process.argv.indexOf("--table") !== -1) {
    logNotice(`chr${CHROMASOME} has ${originalCount} markers, trimming to ${MARKERS_PER_CHROMASOME}, minFreq=${_.min(acceptableMarkers, m => m.minFreq).minFreq} maxFreq=${_.max(acceptableMarkers, m => m.maxFreq).maxFreq}`);
    logNotice(`SLUG\tCHROMASOME\tAVG\tMIN\tMAX\tPOPS...`);
    for (let marker of markersToDisplay) {
      logNotice(`${marker.slug}\t${CHROMASOME}\t${Math.round(marker.avgFreq * 1000) / 1000}\t${marker.minFreq}\t${marker.maxFreq}\t${marker.populationFreqs.join("\t")}`);
    }
  }
  
  if (process.argv.indexOf("--probs") !== -1) {
    let failProbability = 1;
    for (let marker of markersToDisplay) {
      failProbability *= 1 - marker.avgFreq;
    }
    logNotice(`Chr ${CHROMASOME} match probability ${1 - failProbability}`);
  }

  if (process.argv.indexOf("--nojson") === -1) {
    let markers: { [rsId: string]: string } = {};
    for (let marker of markersToDisplay) {
      markers[marker.rsid] = marker.allele;
    }
    SCRIPT_RESULTS.push({
      chr: CHROMASOME,
      markers: markers
    });
  }
}

function printResults() {
  if (SCRIPT_RESULTS.length) {
    console.log(JSON.stringify(SCRIPT_RESULTS, null, "\t"));
  }
}

export interface ChromosomeMarkerList {
  /**
   * Chromosome name: "1" to "22" or "X"
   */
  chr: string;

  /**
   * Map of SNP rsId to allele, e.g. {"rs123": "A"}
   */
  markers: { [rsId: string]: string };
}

interface Marker {
  slug: string;
  chr: string;
  rsid: string;
  allele: string;
  populations: number;
  populationFreqs: number[];
  minFreq: number;
  avgFreq: number;
  maxFreq: number;
}

enum HapmapFields {
  RSID = 0,
  CHROM,
  POS,
  STRAND,
  BUILD,
  CENTER,
  PROTLSID,
  ASSAYLSID,
  PANELLSID,
  QC_CODE,
  REFALLELE,
  REFALLELE_FREQ,
  REFALLELE_COUNT,
  OTHERALLELE,
  OTHERALLELE_FREQ,
  OTHERALLELE_COUNT,
  TOTALCOUNT
}

function fatalError(message: any) {
  logNotice("Fatal error: " + message);
  process.exit();
}

function logNotice(message: any) {
  console.error(message);
}