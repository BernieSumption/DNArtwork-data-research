/// <reference path="../typings/index.d.ts" />

/**
 * Parse a collection of gzip'd hapmap allele frequency files and output a
 * list of alleles that serve as relatively rare markers, i.e. they are within
 * a defined range of frequencies in every popultion.
 */

const MIN_MARKER_FREQ = 0.02;
const IDEAL_MARKER_FREQ = 0.05;
const MAX_MARKERS_PER_CHROMASOME = 1000;
const POPULATIONS = ["ASW", "CEU", "CHB", "CHD", "GIH", "JPT", "LWK", "MEX", "MKK", "TSI", "YRI"];
const CHROMASOMES = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"];
const FREQS_FILE_PATTERN = "./__tmp_hapmap-freqs/allele_freqs_chr{CHR}_{POP}_r28_nr.b36_fwd.txt.gz";


import _ = require("underscore");
import fs = require("fs");
import zlib = require("zlib");
import readline = require("readline");

let analyser: Promise<any>[] = [];

let CHROMASOME: string;
let MARKERS: { [slug: string]: Marker };

let chrIndex = 0;

recursiveAnalyseChromasomes()
  .then(() => logNotice("DONE!"));

function recursiveAnalyseChromasomes(): Promise<{}> {
  CHROMASOME = CHROMASOMES[chrIndex++];
  if (CHROMASOME) {
    logNotice(`Starting analysis of chromasome ${CHROMASOME}`);
    MARKERS = {};
    return chromasomeAnalyser(CHROMASOME)
      .then(printResults)
      .then(recursiveAnalyseChromasomes);
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
          processMarker(cells[HapmapFields.RSID], cells[HapmapFields.REFALLELE], parseFloat(cells[HapmapFields.REFALLELE_FREQ]));
          processMarker(cells[HapmapFields.RSID], cells[HapmapFields.OTHERALLELE], parseFloat(cells[HapmapFields.OTHERALLELE_FREQ]));
        }
      })
      .on("error", reject)
      .on("close", () => {
        logNotice(`   Finished processing ${chromasome}>${population}`);
        resolve();
      });
  });

  function processMarker(rsid: string, allele: string, freq: number) {
    let slug = `${rsid}/${allele}`;
    let marker = MARKERS[slug];
    if (!marker) {
      marker = MARKERS[slug] = {
        slug: slug,
        populations: 0,
        minFreq: 1,
        maxFreq: 0
      };
    }
    ++marker.populations;
    if (marker.minFreq > freq) {
      marker.minFreq = freq;
    }
    if (marker.maxFreq < freq) {
      marker.maxFreq = freq;
    }
  }
}

function printResults() {
  let populationCount = POPULATIONS.length;
  let acceptableMarkers: Marker[] = [];
  for (let slug in MARKERS) {
    let marker = MARKERS[slug];
    if (marker.minFreq > MIN_MARKER_FREQ && marker.populations === populationCount) {
      acceptableMarkers.push(marker);
    }
  }
  acceptableMarkers = _.sortBy(acceptableMarkers, m => Math.abs(IDEAL_MARKER_FREQ - m.maxFreq) + Math.abs(IDEAL_MARKER_FREQ - m.minFreq));
  acceptableMarkers.splice(MAX_MARKERS_PER_CHROMASOME);
  logNotice(`chr${CHROMASOME} has ${acceptableMarkers.length} markers minFreq=${acceptableMarkers[0].minFreq} maxFreq=${acceptableMarkers[acceptableMarkers.length - 1].maxFreq}`);
  for (let marker of acceptableMarkers) {
    console.log(marker.slug);
  }
}

interface Marker {
  slug: string;
  populations: number;
  minFreq: number;
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