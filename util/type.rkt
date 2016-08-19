#lang typed/racket

(provide Dna
         DnaAtom)

(define-type DnaAtom (U 'A 'T 'C 'G))
(define-type Dna (Vectorof DnaAtom))
