# nextflow-base

[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](http://nextflow.io)
[![CircleCI](https://circleci.com/gh/codingene/nextflow-base.svg?style=svg)](https://circleci.com/gh/codingene/nextflow-base)
[![Build Status](https://img.shields.io/travis/codingene/nextflow-base.svg?logo=travis)](https://travis-ci.org/codingene/nextflow-base)
[![Build Status](https://github.com/codingene/nextflow-base/workflows/nextflow-base%20CI/badge.svg)](https://github.com/codingene/nextflow-base/actions)
[![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/codingene/nextflow-base)](https://hub.docker.com/r/codingene/nextflow-base)


## About 

A basic boilerplate structure for [Nextflow](https://nextflow.io) based pipelines with continues integration test case added. This is mostly inspired by [nf-core](https://nf-co.re/) code and structure style.

Currently on this pipeline three basic steps on any Sequence based analysis (starts from fastq files) present -

* Quality Check (using fastqc)
* Filtering (Using fastp) 
* Sequence Read Quantification (Using kallisto)

This can be used as a base to add other `process`.

## Clone and Modify 

```
git clone https://github.com/codingene/nextflow-base.git
```

## Test Run

Without cloning if just want to test run it.

```bash
nextflow run codingene/nextflow-base -profile test,docker
```

This will directly pull the repo from GitHub and execute. 

Check the [Documentation](docs/index.md) for more.
