#! /bin/bash

getorf -sequence $1 -outseq intermediate.orf
patmatmotifs -sequence intermediate.orf -outfile $2
rm intermediate.orf
