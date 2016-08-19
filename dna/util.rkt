#lang typed/racket

(provide reverse-complement
         hamming-distance)

(: reverse-complement (-> Dna Dna))
(define (reverse-complement dna)
  (list->vector
   (map (lambda ([atom : DnaAtom]) (case atom [(A) 'T] [(T) 'A] [(C) 'G] [(G) 'C]))
        (reverse (vector->list dna)))))

(: hamming-distance (-> Dna Dna Integer))
(define (hamming-distance dna1 dna2)
  (vector-count (lambda (a b) (not (eq? a b))) dna1 dna2))
