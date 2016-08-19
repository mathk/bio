#lang typed/racket

(require "../util/type.rkt"
         "../util/search.rkt")


(define (slice-at [list : Dna]
               [fromPos : Integer]
               [toPos : Integer]) : Dna
  (vector-take  (vector-drop list fromPos) (+ (- toPos fromPos) 1)))

(define (slice [dna : Dna]
               [fromPos : Integer]
               [length : Integer]) : Dna
  (slice-at dna fromPos (- (+ fromPos length) 1)))

(: pattern-match-indexes (->* (Dna Dna) (Integer) (Listof Integer)))
(define (pattern-match-indexes pattern dna [pos 0])
  (let ([next-pos (kmp-search pattern dna pos)])
    (if (= next-pos -1)
        '()
        (cons next-pos (pattern-match-indexes pattern dna (+ next-pos 1))))))

(: pattern-count (-> Dna Dna Integer))
(define (pattern-count [pattern : Dna] [dna : Dna])
  (length (pattern-match-indexes pattern dna)))


