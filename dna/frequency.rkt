#lang typed/racket

(provide frequent-word)

(require "../util/type.rkt"
         "../util/database.rkt"
         "util.rkt"
         )

(define-type Hash-Slab (Pairof Integer Integer))
(define-type Frequent-Word (Pairof Dna Integer))


(: frequent-word (-> Dna Positive-Integer (Listof Frequent-Word)))
(define (frequent-word dna length)
  (let ([vect ((inst make-vector (U Zero Hash-Slab)) (expt 4 length))]
        [dna-length (vector-length dna)]
        [vector-filter-t (inst vector-filter (U Zero Hash-Slab) Hash-Slab)]
        [sort-t (inst sort (Pairof Integer Integer) Integer)]
        [vector-set-t! (inst vector-set! (U Zero Hash-Slab))])
    (define (frequent-word-iter [pos : Integer]) : (Listof Frequent-Word)
      (if (> pos (- dna-length length))
          (let* ([sorted-list (sort-t (vector->list (vector-filter-t pair? vect)) > #:key cdr)]
                 [max-count (cdr (first sorted-list))])
            (map (lambda  ([slab : Hash-Slab]) : Frequent-Word
                   (cons (index-to-dna (car slab) length) (cdr slab)))
                 (takef sorted-list (lambda ([k-mer-spec : Hash-Slab]) (>= (cdr k-mer-spec) (- max-count 1))))))
          (let* ([index (dna-to-index dna pos length)]
                 [count (vector-ref vect index)])
            (if (pair? count)
                (vector-set-t! vect index (cons index (+ 1 (cdr count))))
                (vector-set-t! vect index (cons index 1)))
            (frequent-word-iter (+ 1 pos)))))
    (frequent-word-iter 0)))

(define (neighbors [pattern : Dna]
                   [distance : Integer]) : (Listof Dna)
  (let ([pattern-length (vector-length pattern)]
        [current-neighbors ((inst make-hash Integer (Listof Dna)) (list (cons 0 (cons pattern '()))))]
        [alter-nucleotide ((inst make-hash DnaAtom (List DnaAtom DnaAtom DnaAtom)) '( (A . (T C G)) (T . (A C G)) (C . (A T G)) (G . (A T C))))])
    (let neighbors-iter ([pos : Integer 0])
      (if (>= pos pattern-length)
        (set->list (list->set (apply append (hash-values current-neighbors))))
        (begin
          (for* ([key (hash-keys current-neighbors)])
            (let* ([new-key (+ 1 key)]
                 [value (hash-ref current-neighbors key)]
                 [new-value (hash-ref current-neighbors new-key (lambda [] : (Listof Dna) '()))])
            (when (< key distance)
              (hash-set! current-neighbors new-key 
                         (append (apply append (map (lambda ([previous-neighbor : Dna]) : (Listof Dna)
                                                      (map (lambda ([nucleotide : DnaAtom]) : Dna
                                                             (let ([new (vector-copy previous-neighbor)])
                                                               (vector-set! new pos nucleotide)
                                                               new)) (hash-ref alter-nucleotide (vector-ref pattern pos))
                                                      )) value))
                                                    new-value)))))
            (neighbors-iter (+ 1 pos)))))))

(format-dnas (neighbors (string->dna "TAAGTGTTGA") 2))



