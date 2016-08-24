#lang typed/racket

(provide kmp-search
         kmp-table)

(: kmp-table (All (A) (-> (Vectorof A) (Vectorof Integer))))
(define (kmp-table pattern)
  (let* ([t (make-vector (vector-length pattern) 0)])
    (define (kmp-table-iter
           [pos : Integer]
           [cnd : Integer]) : (Vectorof Integer)
      (vector-set! t 0 -1)
      (vector-set! t 1 0)
      (cond
        [(>= pos (vector-length pattern)) t]
        [(eq? (vector-ref pattern (- pos 1)) (vector-ref pattern cnd))
         (vector-set! t pos (+ cnd 1))
         (kmp-table-iter (+ pos 1) (+ cnd 1))]
        [(> cnd 0) (kmp-table-iter pos (vector-ref t cnd))]
        [else
         (vector-set! t pos 0)
         (kmp-table-iter (+ pos 1) cnd)]))
    (kmp-table-iter 2 0)))

; Perform a searh in a vector using the
; Knut-Morris-Prat algorithm.
(: kmp-search (All (A) (->* ((Vectorof A) (Vectorof A)) (Integer) Integer)))
(define #:forall (A)
        (kmp-search [pattern : (Vectorof A)]
                    [dna : (Vectorof A)]
                    [next-offset : Integer 0]) : Integer
    (let ([dna-length (vector-length dna)]
          [pattern-length (vector-length pattern)]
          [jump (kmp-table pattern)])
      (define (kmp-search-iter [offset-dna : Integer] 
                               [offset-pattern : Integer]) : Integer
        (cond
          [(>= (+ offset-dna offset-pattern) dna-length) -1]
          [(eq? (vector-ref pattern offset-pattern) (vector-ref dna (+ offset-dna offset-pattern)))
           (if (= offset-pattern (- pattern-length 1))
             offset-dna
             (kmp-search-iter offset-dna (+ offset-pattern 1)))]
          [(> (vector-ref jump offset-pattern) -1)
           (kmp-search-iter
             (- (+ offset-dna offset-pattern) (vector-ref jump offset-pattern))
             (vector-ref jump offset-pattern))]
          [else (kmp-search-iter (+ offset-dna 1) 0)]))
      (kmp-search-iter next-offset 0)))

