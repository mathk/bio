#lang typed/racket

(provide reverse-complement
         hamming-distance)

(require "../util/type.rkt")

(: reverse-complement (-> Dna Dna))
(define (reverse-complement dna)
  (list->vector
   (map (lambda ([atom : DnaAtom]) (case atom [(A) 'T] [(T) 'A] [(C) 'G] [(G) 'C]))
        (reverse (vector->list dna)))))

(: hamming-distance (-> Dna Dna Integer))
(define (hamming-distance dna1 dna2)
  (vector-count (lambda (a b) (not (eq? a b))) dna1 dna2))

(: hamming-jump (-> Dna (Vectorof Integer)))
(define (hamming-jump pattern)
  (let iter ([out (make-vector (vector-length pattern) 0)]
             [pos 0]
             [hamming 0])
    (vector-set! out pos hamming)
    (if (>=  (+ 1 pos) (vector-length pattern))
      out
      (if (eq? (vector-ref pattern pos) (vector-ref pattern (+ 1 pos)))
          (iter out (+ 1 pos) hamming)
          (iter out (+ 1 pos) (+ 1 hamming))))))
(: min-hamming-distance (-> Dna Dna Integer))
(define (min-hamming-distance dna pattern)
  (let ([pattern-length (vector-length pattern)]
        [dna-length (vector-length dna)])
    (let iter ([offset-dna 0]
               [offset-pattern 0]
               [current-distance 0]
               [current-min-distance pattern-length])
      (cond
        [(> (+ offset-dna pattern-length) dna-length) current-min-distance]
        [(>= current-distance current-min-distance) (iter (+ offset-dna 1) 0 0 current-min-distance)]
        [(eq?  (vector-ref pattern offset-pattern) (vector-ref dna (+ offset-dna offset-pattern)))
         (if (= offset-pattern (- pattern-length 1))
             (iter (+ offset-dna 1) 0 0 (min current-min-distance current-distance))
             (iter offset-dna (+ 1 offset-pattern) current-distance current-min-distance))]
        [(= offset-pattern (- pattern-length 1)) (iter  (+ offset-dna 1) 0 0 (min (+ 1 current-min-distance) current-distance))]
        [else (iter offset-dna (+ offset-pattern 1) (+ 1 current-distance) current-min-distance)]))))

#|
(: min-hamming-distance (-> Dna Dna Integer))
(define (min-hamming-distance dna pattern)
  (let ([pattern-length (vector-length pattern)]
        [dna-length (vector-length dna)]
        [jump (hamming-jump pattern)])
    (let iter ([offset-dna 0]
               [offset-pattern 0]
               [missmatch 0]
               [current-distance 0]
               [current-min-distance pattern-length])
      (printf "off-dna ~a, off-pat ~a, missmatch ~a, current ~a, min ~a~n" offset-dna offset-pattern missmatch current-distance current-min-distance)
      (cond
        [(> (+ offset-dna pattern-length) dna-length) current-min-distance]
        [(>= current-distance current-min-distance) (iter (+ offset-dna 1) missmatch  0 (vector-ref jump missmatch) current-min-distance)]
        [(eq?  (vector-ref pattern offset-pattern) (vector-ref dna (+ offset-dna offset-pattern)))
         (if (= offset-pattern (- pattern-length 1))
             (iter (+ offset-dna 1) missmatch 0 (vector-ref jump missmatch) current-min-distance)
             (iter offset-dna (+ 1 offset-pattern) missmatch current-distance current-min-distance))]
        [(= offset-pattern (- pattern-length 1)) (iter  (+ offset-dna 1) missmatch 0 (vector-ref jump missmatch) (min current-min-distance current-distance))]
        [(= 0 missmatch) (iter offset-dna (+ offset-pattern 1) offset-pattern (+ 1 current-distance) current-min-distance)]
        [else (iter offset-dna (+ offset-pattern 1) missmatch (+ 1 current-distance) current-min-distance)]))))
|#
;; Adjust next mismatch

