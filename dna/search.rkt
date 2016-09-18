#lang typed/racket

(require "../util/type.rkt"
         "../util/search.rkt"
         "../util/database.rkt"
         "util.rkt"
         math)

(require/typed racket/function
               [(curryr curryr2) (All (a b c) (-> (-> a b c) b (-> a c))) ]
               [(curryr curryr3) (All (a b c d) (-> (-> a b c d) b c (-> a d)))])

(current-logger (make-logger 'dna-search (current-logger)))

(define-struct/exec profile-unit ([a : Exact-Rational] [t : Exact-Rational] [c : Exact-Rational] [g : Exact-Rational])
                    [(lambda ([pu : Profile-Unit] [nucleide : DnaAtom])
                       (case nucleide
                         [(A) (profile-unit-a pu)]
                         [(T) (profile-unit-t pu)]
                         [(C) (profile-unit-c pu)]
                         [(G) (profile-unit-g pu)])) : (-> Profile-Unit DnaAtom Exact-Rational)]
  #:type-name Profile-Unit)

(define-struct frequency-unit ([a : Integer] [t : Integer] [c : Integer] [g : Integer])
  #:mutable
  #:transparent
  #:type-name Frequency-Unit)

(define (most-frequent-nucleide [pu : Frequency-Unit]) : DnaAtom
  (let*-values ([([first : DnaAtom] [max1 : Integer]) (if (>= (frequency-unit-a pu) (frequency-unit-t pu))
                                                          (values 'A (frequency-unit-a pu))
                                                          (values 'T (frequency-unit-t pu)))]
                [([second : DnaAtom] [max2 : Integer]) (if (>= max1 (frequency-unit-c pu))
                                                           (values first max1)
                                                           (values 'C (frequency-unit-c pu)))])
    (if (>= max2 (frequency-unit-g pu))
        second
        'G)))

(define (frequency-add! [fu : Frequency-Unit] [nucleide : DnaAtom]) : Void
  (case nucleide
    [(A) (set-frequency-unit-a! fu (+ 1 (frequency-unit-a fu)))]
    [(T) (set-frequency-unit-t! fu (+ 1 (frequency-unit-t fu)))]
    [(C) (set-frequency-unit-c! fu (+ 1 (frequency-unit-c fu)))]
    [(G) (set-frequency-unit-g! fu (+ 1 (frequency-unit-g fu)))]))

(define (frequency-sum [fu : Frequency-Unit]) : Integer
  (+ (frequency-unit-a fu) (frequency-unit-t fu) (frequency-unit-c fu) (frequency-unit-g fu)))

(define (frequency-unit->profile-unit [fu : Frequency-Unit]) : Profile-Unit
  (let ([sum (frequency-sum fu)])
    (profile-unit (/ (frequency-unit-a fu) sum)
                  (/ (frequency-unit-t fu) sum)
                  (/ (frequency-unit-c fu) sum)
                  (/ (frequency-unit-g fu) sum))))

(define-type Profile (Listof Profile-Unit))

(define-type Motifs (Listof Dna))

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

(: median (-> (Listof Dna) Integer Dna))
(define (median dnas length)
  (caar ((inst sort (Pairof Dna Integer) Integer) 
         (map (compose (lambda ([dna : Dna]) (cons dna (sum-min-hamming-distance dnas dna))) (curryr2 index-to-dna length)) 
              (range (expt 4 length)))
               <
              #:key cdr)))

(define (most-probable-k-mer [dna : Dna] [profile : Profile]) :  Dna
  (let ([profile-length (length profile)]
        [dna-length (vector-length dna)])
    (let iter ([pos : Integer 0]
               [best-k-mer-pos : Integer 0]
               [best-match : Exact-Rational 0])
      (log-info "most-probable-k-mer pos: ~a, best match ~a" pos best-match)
      (if (> (+ pos profile-length) dna-length)
          (vector-copy dna best-k-mer-pos (+ best-k-mer-pos profile-length))
          (let ([match (for/product : Exact-Rational ([dna-pos (range pos (+ pos profile-length))]
                                     [unit profile])
                         (unit (vector-ref dna dna-pos)))])
            ;;(printf "current match ~a, best ~a, current pos ~a~n" (fl match) (fl best-match) pos)
            (if (> match best-match)
                (iter (+ 1 pos) pos match)
                (iter (+ 1 pos) best-k-mer-pos best-match)))))))

(define prof (list
    (profile-unit 364/1000 333/1000 182/1000 121/1000)
    (profile-unit 333/1000 182/1000 182/1000 303/1000)
    (profile-unit 303/1000 303/1000 212/1000 182/1000)
    (profile-unit 212/1000 212/1000 303/1000 273/1000)
    (profile-unit 121/1000 364/1000 182/1000 333/1000)
    (profile-unit 242/1000 152/1000 303/1000 303/1000)))

(define my-dna (string->dna "TGCCCGAGCTATCTTATGCGCATCGCATGCGGACCCTTCCCTAGGCTTGTCGCAAGCCATTATCCTGGGCGCTAGTTGCGCGAGTATTGTCAGACCTGATGACGCTGTAAGCTAGCGTGTTCAGCGGCGCGCAATGAGCGGTTTAGATCACAGAATCCTTTGGCGTATTCCTATCCGTTACATCACCTTCCTCACCCCTA"))

;(most-probable-k-mer my-dna prof)

(define (build-profile [motifs : Motifs] [fu-unit : (-> Frequency-Unit) (thunk (frequency-unit 0 0 0 0))]) : Profile
  (for/list ([pos (range (vector-length (first motifs)))])
    (let ([fu (fu-unit)])
      ((inst for-each Dna) (compose ((inst curry Frequency-Unit Void DnaAtom) frequency-add! fu)
                                    ((inst curryr2 Dna Integer DnaAtom) vector-ref pos)) motifs)
      (frequency-unit->profile-unit fu))))

(define (motif-consensus [dnas : Motifs]) : Dna
  (let* ([dna-length : Integer (vector-length (first dnas))]
         [frequency : (Vectorof Frequency-Unit) (make-vector dna-length (frequency-unit 0 0 0 0))])
    (for ([pos (range dna-length)])
      (vector-set! frequency pos  (frequency-unit 0 0 0 0))
      (map ((inst compose Dna DnaAtom Void) (curry frequency-add! (vector-ref frequency pos)) ((inst curryr2 Dna Integer DnaAtom) vector-ref pos)) dnas))
    (log-info "motif-consensus: freq-vect ~a" frequency)
    (vector-map most-frequent-nucleide frequency)))

(define (motif-score [dnas : Motifs]) : Integer
  (sum-hamming-distance dnas (motif-consensus dnas)))

(define (k-mers-length [dnas : (Listof Dna)]) : Integer
  (vector-length (first dnas)))

(define (greedy-search [dnas : (Listof Dna)] [length : Integer]) : (Values Dna Motifs)
  (let ([dna-length : Integer (k-mers-length dnas)])
    (let iter ([best-motif : Motifs (map ((inst curryr3 Dna Integer Integer Dna) vector-copy 0 length) dnas)]
               [pos : Integer 0])
      (log-info "greedy-search best-motif: ~a~ngreedy-search pos: ~a" best-motif pos)
      (if (< dna-length (+ pos length ))
          (values (motif-consensus best-motif) best-motif)
          (let* ([motif : Motifs (list (sub-dna-at (first dnas) pos length))]
                 [profile (build-profile motif (thunk (frequency-unit 1 1 1 1)))])
            (log-info "greedy-search first motif ~a" motif)
            (map (lambda ([dna : Dna])
                   (set! motif (append motif (list (most-probable-k-mer dna profile))))
                   (set! profile (build-profile motif (thunk (frequency-unit 1 1 1 1))))) (cdr dnas))
            (log-info "greedy-search new motifs: ~a" motif)
            (if (< (motif-score motif) (motif-score best-motif))
                (iter motif (+ 1 pos))
                (iter best-motif (+ 1 pos))))))))

(: optimize-search (->* (Motifs (Listof Dna)) (Integer) Motifs ))
(define (optimize-search [motif : Motifs] [dnas : (Listof Dna)] [score : Integer (motif-score motif)]) : Motifs
  (let* ([profil (build-profile motif (thunk (frequency-unit 1 1 1 1)))]
         [next-motifs (map (curryr2 most-probable-k-mer profil) dnas)]
         [next-score (motif-score next-motifs)])
    (if (< next-score score)
        (optimize-search next-motifs dnas next-score)
        motif)))
  
(define (random-search [dnas : (Listof Dna)] [length : Integer]) : Motifs
  (let* ([dna-length : Integer (k-mers-length dnas)]
        [gen : Pseudo-Random-Generator (make-pseudo-random-generator)]
        [gen-motif : (-> Dna Dna) (curryr3 sub-dna-at (random (- dna-length length) gen) length)])
    (let iter-time ([count 10000]
                    [best-motif (map (curryr3 sub-dna-at (random (- dna-length length) gen) length) dnas)])
      (let* ([motif (optimize-search  (map gen-motif dnas) dnas)]
             [keep-motif (if (< (motif-score motif) (motif-score best-motif)) motif best-motif)])
        (if (<= count 0)
            keep-motif
            (iter-time (- count 1) keep-motif))))))
  

;; (motif-consensus '(#(G C G T) #(A T G T) #(A C G A)))
;; (require "../util/database.rkt")
;;(define dnas (strings->dnas (list "AAATTGACGCAT" "GACGACCACGTT" "CGTCAGCGCCTG" "GCTGAGCACCGG" "AGTACGGGACAG")))
;; (median dnas 3)
;; (map (lambda ([dna : Dna]) (min-hamming-distance dna #(C A A))) dnas)

;; (define dnas (strings->dnas (list "AGTTGTGTGAGCGCAGCCACGCTCAAATAGCTAGTGGAAGAA" "AGGACCTCGGGTGCATGCCCTCGAATAGTGGGGGTAAACGTA" "CTGTCGAACTCAACAGGTATCCAGGTACTTGCCCTACTAGTG"
;;                                   "ATACCGAGTCGCACTGTATGGTGAAATCAGGTAGTGGCACAA" "CCCACTGTGCAATATATGTTTAAGATGATCTTAGTGAGTACC" "AGGAATTGACCACTAGTGCATCCATACTTACACGTGAGATAT"
;;                                   "GTGGTCGTCAGGCCTCAGTTATGTATGAACTTAGTGGTGCTT" "CTGCTTAAGCAGGTAGACGCTTCCCGCAATCTCAGCATAGTG" "ACCAATGCTACACGTTGCACAGCTATAGTGATCTGACACCTT"
;;                                   "GAATGTATGCCAGAAACGTCCCTAAGTAGAGATTTTGTAGTG")))

;;  (greedy-search dnas 5)
(define dnas (strings->dnas (list "GGCGTTCAGGCA" "AAGAATCAGTCA" "CAAGGAGTTCGC" "CACGTCAATCAC" "CAATAATATTCG")))
;(greedy-search dnas 3)

(define dnas-sample (strings->dnas (list
                                    "AGGCGGCACATCATTATCGATAACGATTCGCCGCATTGCC"
"ATCCGTCATCGAATAACTGACACCTGCTCTGGCACCGCTC"
"AAGCGTCGGCGGTATAGCCAGATAGTGCCAATAATTTCCT"
"AGTCGGTGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTG"
"AACCGGACGGCAACTACGGTTACAACGCAGCAAGAATATT"
"AGGCGTCTGTTGTTGCTAACACCGTTAAGCGACGGCAACT"
"AAGCGGCCAACGTAGGCGCGGCTTGGCATCTCGGTGTGTG"
"AATTGAAAGGCGCATCTTACTCTTTTCGCTTTCAAAAAAA"
)))

(define dnas-sample2 (strings->dnas "CTCACTCACGCGCCGTAGGTTAAGATGCACTGGTGGCTCTATATCGATGTTGCTTCTTTCAGCCCCGAAGGGCCTGACTCTCCTTGTCGCCATGATTTTAACTTTTATTGATGACGCCGAAACAAGATCCAGTCCAGTTTTGCAAGCGACAGCAAA"
"TCTCACTCCAATGACACCTATGTCGATGATGCTTGCACGTTTCAGAGTATGAGGCCGTCCGGATTGAGCCGAGATATCACCCAGTAAATTGAGATCGTGCGGTCTTGCGGACATAACACTGAAGAGGTTCTATCTGCACGCCACGTACTATCGCTA"
"AGATATTGAACAAGCGGGTGGGTCAGCCGGGACGTTCTCCATATTACCATGCAATCAGTCGCCAACATAACGCTGACCACGTACGATGATACTTGCAAAATCCATCCACACGCCGGGGATCGGGTGGACGACATATAGTTAGGACGATCGTGCTCC"
"GGGGAACTTCGGCCCTTTGCCGTTCTCATCACGCACGGTGTTACTTACAGGAACGACATCATCCGATGCCGCCCGTATGCAGTGTTGCCCGCTCCACTTCGCTGGCCATCCACACGACCCTATGTGTAACGTCTTGAGATGAGATACGGATGACAT"
"GCTGATGCTTACCGCGTGGGCGCACCCTCCTCGCCCCCTTTTCGCACACTGCGCGAGAGTTCCTCGATTTGAATCCGTAAGTACAGCGCGGCCTGGCTCTCTCAGTGGACTACACCGACTGAATGGATTAACCAATGGCGTAACGGTGTAGAAAGA"
"GTTGATACTTTCATCGCATCGTGCCGGGTGGAAGCTAGGTTCCGTTAACCCCTCACCCCTACTCAAGTTCTTCATCGATTCATACCTTTGTACGAACTTTGCTGTAAGTCAGGTGTGCAAGGGTGGGCAATATCAATATGTATTAAACCACAGCCC"
"CAGAAGGTGATGGTTGCTGCTTACACAGCGACGTCGTCCCCCCGCTGCTTGGGCGGGAGTAGCCTATGCGGTCGAGAGGTCGAGTGGTGGCATTAGGACAAAGGGTGCCCGAAGTGACTCTCATGCACCCCTCTGGCAGATGGAACAAGATAAAGA"
"CAAAATAATTATGTTGGTCCTTACTTTCAGTCGGCAAACCTACCGTGCTATTCTCCCCTAGTTGCGAAGAGCTCTCTATTTAGGCCGTAACGTACACCCCTGGGATATCAACCAACTCAAAGGTCTAACCTCATTGGGGCCAAACAGATGGGACCG"
"TCGCCAAACCAAGAATATTCTCTGACTCTGAGTCCTACACGCCGGACCAGAAAGGGCACTTAAAGTCGAACAGATGGTTCTTTCGTAGCCTTTTCAGTATGCATCTACTCAGGTGAGGGCTACAATTATACGTTTTCGGTTAAGGTCGGCTTCAAT"
"AGTCCACACCATGAAGTCGCCGTCCCCCCCACTGATTGCGGATAATGTTGCTGCAAAATCGTTGGTTCTTGCCTAATATCGGGTATATCAACTAGGCTGCGCAGATCCATTTAAACTGACTAGGTGGACTACTTATCCATACCGCGCAAACGATCG"
"CGGTTCGGTAAAATTTCCAAAGACCCCTAGGCAAGAGTTGCTCCTTCCGAGTGCTAGATCAGGCGGCAATTGGAGACAGTGGAACAGTTAGACGAGGCCTGGAGCATGGCCGAGGCGGAACCACTTCCATCAATCGGAGAATACAGACACACGGAA"
"AAGAAGTATGAGACTACGTTTACATTCTAATGTGACATGGGATACGCCGGAAGTAAATTGGAACTGGAGTTTAGCCGCGCCGGAACTGTGGGTTGCGATGGTCCTTACTGACGGCTTACACCATGATCCGATGTGTCTCAGCAGAAACTGCTCCAA"
"ACACACGATTCAACACTGCCGTTCCAGAAGCATAGAGGTTCAGCCGTCACCCAGACAAGCTTTGGATCTCAGATCATGCTGCCGGTTGATTCTTCCGGGCGTAATAGAGCCGAAATACACTCAAGGATTAACTCTATTTCAAGGGCTCTCTCGTTT"
"AGAGGCATGGCCGGTATCCTGTCGTACATGAAACTAGTAGGGTTTGCCTCAAGTAATGGAATACTTGGTGTTCTCGGTTAAGGGTGACGCGAGGGGCGACCAAACGTGGGTGCTTCTTACCCCCGACCCGCTCCGTGTAGTGCGCCCGGTTGGTTA"
"AAGCCGGCATCGTATTTATCTAAAGTTGTTGCTTGCCAGTTAGATACCGTACGCCTAGGGAGGCGGTACGCATGTACAACATTCGATTGCCACCACCAAACGCTAGCTGTGACGCTTTCATATGTCGCTCCACACGGCCGGATAGACGCGATACAT"
"TGCACTGCTACATTGTCAAATGAATAGGCAATTGTTCACCTCTTGCTATAACTTCGGCATTGTCGATTCCAATAGCAGCTCAATCACGGGTCAAACAGGAGGAGCGCCGATAAAGGCTGAACCACCACGAAGGATGATACTTGCCCGTCGACAAAA"
"ATGATATTGATTCTCTATAGCAGTCACCGGAGTCGGTAGCATCGACGGGGTGCTCCTTACGGACTCACCACGGAACGGGGTATTACCCCAGTTCTTGAGCGGGTAACCGATAGTCTATCAAGCCAAGAACAAAGCACCAAAACTGAAAACCAAGAC"
"AGTAGTACGCTGTCCGAGTGTTAAATCATTGTGAGTACTGGTCTTTTGAGGCGAGGTCGCCCTATCCCCTAATTTCATAAGAGTTCCCCCGGAATCGATGCTACTTACAGAAAGCGGGACTTCAATCCGCATTTACCGTTACGCACAGATGCGCGC"
"CCTCTTAAAGTAAACCTAGCGGAGATACGCATTGATCCAAGAGAGTCGTGGTCACGATAAGACCGCAGGTCTACAATCGACAAAATTACTTTCTATGCTACACCGGATGATGCTACTTACGCGCACGTTTTTCTTCACCATAGCATTCAGCTGGAC"
"GGCATGCCCACCCTCATAGCCATTATTGCTCGTCTTGTTTCCAGATTCTGCATCTTTTGTTAGTATTGTTCGGTTGATGCTTACATATGGACAACCTTGCTTATTGGAATGGTTATACGTGGAGTGGCGACGTCGGAGTCAGGGCAATTTGTTTAG"
"TCTTGCGCGCGACTCAATACACGAGTTAAGCGACCATAGGCACTACGCCTTCAAGGGGAGACACGTAAACCTGAGAGTTCGCGGGGGAAGAACGAAGGTCACGGGATTCCGTTCCGAGCCGGTGGTTCTTGCTTGTGAAAGTACGACAGCTAAGGG"
"AAGTCACCCCGAGGCCATGTGGGCGGATCTGAACTTAACCGTTATGTGTCCTACACGGCCAATTGTGCGTATGCTGCTACTTGCCTACATTAGAGCCAGTATTACTGACCAGAGGTGAAAAATAATAGTTCTGCTATACACTTGAAAGAAAGACAC"
"GAGTGGATTTAAGTTACCACAGAGTGTGTAGCTATACGCTGTGTAATGTCCCACGTATTGGAGTTGAATTCTGATCCATTGCTAAGAACATCCGTGACTAACTCCATCAGACCTGGGTTAGACAAGTATTCAGTTGCTCCTTACGTTTTTGCAACG"
"CATTGTTTCTAGTCTACCCATTTTTTGAGTGCTGCGTCATGATACCAACCCGGAACGTTATCGGACCCGGACCAGACAGTTACAATTGAGCCTCAGAGTCTCAATATATTGATACGTAGTGTCACACTGAGGGGTGTTGCTTCCTTTCCCAATACG"
"TCGTCTCCCCCCGGACACTATGCGCTCTTAGCCATAACGCAGCTCATTATATCTTCGTGCGGTGTTGCTTTCCAACGAAAATTACAAGAGGGCTGCTCTATGAAAAGTCGGAGTGATATTAAATCGACTGAATGTTTCCTGCCGTTCACAGAAAAA"
))

;(define-values (consensus-ex motif-ex) (greedy-search dnas-sample 5))

;(format-dnas motif-ex)
;(printf "~n")
#|
(let-values ([(consensus2-ex motif2-ex) (greedy-search dnas-sample2 12)])
  (format-dnas motif2-ex)
  (motif-score motif2-ex))

 (let ([motif (random-search dnas-sample2 12)])
  (format-dnas motif)
  (motif-score motif)) |#

(define dnas-ex (strings->dnas "CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA"
"GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG"
"TAGTACCGAGACCGAAAGAAGTATACAGGCGT"
"TAGATCAAGTTTCAGGTGCACGTCGGTGAACC"
"AATCCACCAGCTCCACGTGCAATGTTGGCCTA"))

(let ([motif (random-search dnas-ex 8)])
  (format-dnas motif)
  (motif-score motif))
  

;(motif-score (strings->dnas "AGGCG" "ATCCG" "AAGCG" "AGTCG" "AACCG" "AGGCG" "AGGCG" "AGGCG"))
;(motif-score (strings->dnas "AGGCG" "TGGCA" "AAGCG" "AGGCA" "CGGCA" "AGGCG" "AGGCG" "AGGCG"))

(define prof1 (build-profile '(#(A G G C G) #(T G G C A) #(A A G C G) #(A G G C A) #(C G G C A) #(A G G C G)) (thunk (frequency-unit 1 1 1 1))))
(most-probable-k-mer (string->dna "AAGCTTCCAACATCGTCTTGGCATCTCGGTGTGTGAGGCG") prof1)

#|
(define dnas-sample (strings->dnas (list "GCAGCATAGGTACACTCTCGGCAATCGGCCAGCCGTAATTGATTCGTCAACGAAAACGCGAGCCAACTTACTAAGCGGCTACGCGATAATAATTCGCAATAATGCTAACGAAAACGCTTTACAAAGTTGCATTGAATACGATGTATCTCTTTTACC"
"GAGCACAATGCCTCAGAGACCCCAACTACGTTACATCATTCCAGAGGGTCCCGTGGCGTATGGATAGCACAATGCTTTACAGCACAGTAATTTCGGGTGCGTAACCCATGGTCCACCATACGCATTAGGTCGTAGGGGTCGAGTTCATTTATCCTG"
"CTCAGCGCGGGAGGCACCCCTGATGTGGGCAGTCTGAAGAAAAGACAACTGATCCTTGAGATGGGTATGATATGATCATTATACCGACGCAATTAAGTCAGCCGATGAGATGGCGTGTGCTAACCCCGGGCAACAAGGTTTCATCATTTGTTAAAT"
"TAATTTATAATGATAACTGGTTGCTAAGTCGCACGTTTAAAGATGAAATTCGAGACGTTTTATATTGTACATGCTTCATGTGCTTGGCAGATGGGTTCGAGATAATATGTCATTCTCCATTGAGTTCACGAGGTTCGCCTCGTAACGAGGTTGCGT"
"TTTTCTCGTCGCAAGAAAAATAAAAAGAAAAATCATTGTAAAATCCTGGCGGAACCGGCTTTGCTCATTAGAACCAAGTTCCTTGATTGGTCCAAGGAGGCGTGCTCGCCTGGACAGTTTTTTAAAGAGGCAGACGTCTAAGGCTTCGTCTTGGTT"
"CAGGGTTACCCCGAAGCAACGCCGCGCAAACCATGACAAGCGCTACCTCTCGATCCCTGCCGGCCGAGATAACTACCCGCCTACACTAAGTTACATGTAAACTTTGGAGGCATTCTCATTCCGCCCATTGGCCGCATTTAGACTACCCAGTATGCC"
"GAAAACCCGTACTTTGGAGGGATTTACACCCAATCGAGCCCTGTGGTTTGACCTTTTCAAACCAAGTTCCGTCAGTATTTATCGGCATACTATTGCTGCGTAGAATAGTATTAGGTGGGGCTTCCGGGACAGGTGCTCCGCGTGTGCAACACTACG"
"ACCATGTTGCGTAGCCATATACCTATAAGGATTCGCCGGATCGTGACCAGTCGCCTGCCTACTTTCGTCAGGTCCATAGGTACGGCTGTCTCTCGAGTATAGCCGCCAGTCCCCGGTCGCACGGACTTATCAGCTCATAGATAACGGCCATAGACT"
"CTGGGTAAGCTCTGGTGGTCGGGAGACCTATGATTGAGCACAGGTGAAGGGGCTCCGGAGCTGCTTACCCTTACCAAGTTGCTTCTCACCACGGGCATGGAACCTGGTGAATCCACAGGTAACTAAACCTAAACCCATATTAGATGTAGATTTCTT"
"TTGTGAGTCTTCTTGTTTTACGGAGGAGGGCGACCGCCCTGCTTTACCGTCCACCACTCCGATTGCGAGTGTGTCCTGTTGGGGGAAGGGTTGTCAAGTCTAATCTCGAACCGCCACCTACATTTGCAGACACGCTACGAAAAGACCACGTTTCAT"
"TCTTTAGTTGTCAACCAACTGCTAACTAAGTTGCGTACCCGTAAAACTTCACGTGACGACGAATTATGTATCGCTGTTTAAATACACCTTTGAGTAGGTGCCTTTGCAAAGACCAGAAGTTATTCCCCGGCAGCTTCATTCGAACGTACGCAACCT"
"GGGCCACTAGAGCAAAGGGGCCTGGGCCGGCGCGGACATAGGTAGGCCGTCGCGTACGACTTCGTGTGGGTTGCCGTGCATCCGGCTACATTAATGGAAGAAGGCTAGACCAGGTTTCCTGCTACTCGAAGGTATAACTTTTAGTAGGATCAATGT"
"GAGATGGTAATAGAACGCTTGCACTCAATGAGGAGTGGCCGTTGTAACACACCTCATACGAATTTAATCGTCGGGCAGCCAAGCACGACGTTACCTCGCGGCCCTGTGACTCCAGCCAGTAACATTTTGCCGTTTGACCTCCTAAGATGGTCATGG"
"GCGTTCAAGGGAGGGAGATAAAATGTTTTTAAGAGTCTGCCAACAACTAAGCCTCACAGATCGAGGTCGGTTACGACGTTACGTCTGTTTCTTTTAACGTCACTGTTAGCGGGCTCCAGAGATTAGGTGGGTTGGCAGAAATCTGCCCGTCAGCGA"
"CGACTGAAATGTGTCACCGTTTGAGGAATTCTGTGTACGAAGTTGCATGTTTCCAGGCACCGCATAGTCTCCGAATATACAAGTATTTTGAGATTATTGGAATCTAGATAGCTTAGAGAATAATCTGCGTCGCTCTCATCCAACGAGTTGGACGTG"
"TCGACTGAACATCATCGGCGTCGTCCATACCCACTCCTCATGTGCTATTTAAACGTTGTACGAGAGTATAAGTCATTGGGGGACGTCACAATTATTAACAAGAAGCACACTAGGTTCCTTGCTTTTATTCGGACAGGTCATACTTGCTTCGCTAGA"
"CTTCTATAGAACCTAATATCAGTTCGGTCTTCAGGGCCGGGAGTTGACGCGGTTCAGCGCACTATCCGAGAATCGGTGAGGAGGACAACGTTGCGTGCCACGTTTATACTGAATTCGATTGTGGATATAAGGGCTTCAGAACGGCGGTTTAATAAA"
"GCATTGGGACTTGTGGTTAGATCCACCATGTTGCGTCAGGGGAGATCCTCCGAAGTCCTCTAACTCTTCAGGGGGGTTTGTCTGGAGCGGCGGGTAATGCTCCAGCAACTGGCGGTAATGCATCCAATGTTAAGGCACAGTCGCGTACCGCACTGT"
"GTCGGCACCGTCGACATACGCAGTCGTGAAACATAGTTCCCTTTATCTTCGGTCCAAGGGCCTGAGTTACTTGGCAATTTATCGCCAAGCGCGCATCAGAAGTGGTAGGGGCTGTCAATTATAAGTGACAATACGGGACATAACACTATGTTGCTT"
"ATTGGCTCTCCCCCGGGACGTTAAGTTCCCTGATGATTCGGAGGCTCGCATCCCAGTCTTCAGAAGCATGGCACGATGTTGCTTGGGACCACTGGGGGGAACGGTCTATCCTTCGCATACGGGGTTAAAGCTCAGGATGTTTACTAGAGCTCTGCG"
"CCAGACTCTAAATGATGCTCAACCTAAACAAACTCTCAGGTACTGGGCGCGTATATCCTGAAATAAGGGCACCTGTTGGATTAACCCGGAAGAGTTCCAGCGCTTTAATGGGAAACGTGAACAATGTTACTTCGTTTGCAAGAAGTAACGTAATCC"
"CACTACCGGTTGGCGAGACTAATCTATAAGTGGGGAACGCCCACGGGTAAGGTCTGTCGTACAAAGTTTCATACATAAGCACCCTTACCGTCCGGCAGAATCAATGTTCCCTCATTGAAACTCTTACCCTAGATGTACTGTTTAAAATCCTGCTTC"
"TATTGGCAGATCTGGATTTAAACTGGTCAATGTTTTACCGCCCTGCACTGGTGCAGGTTGTCCCCCGTAATTTTTCGCCAGCAGACAATGTTACTTGGGTTTCCTTTTCCACAGCACAGGATACTTCTTCTCGAACACAGCTCTTTGATCTGTATG"
"AAATAGCCTCATCCCCGAGCACAGCGGAGCGGGCTTTATCAGATTAGCGTGGTGATAGGCTGAATGCAATAGAAGCAAAGGCTTTTTTCCCACTGAACAATGTTCCGTCTATAACCAGAATCCAGGGTACTAACGCCATGGTTAAATATCCTGTGT"
"AGTATATAATCCTGTGACTTTGTTTTGAAAACAGCGAGTTACCTAGTAGGGCGTACGATCCGTGTACAACCCCACCTTACGCTCCCGGACTGTTCCAACGTGAACAGGTCTAATGGAATATCACTCGGGACAACGAGGTTTCATCTTTTTGTAGCC")))

(define-values (consensus-ex motif-ex) (greedy-search dnas-sample 12))
|#


; (require/typed profile [ profile-thunk (-> (-> Any) Any)])
; (profile-thunk (thunk (random-search dnas-sample2 12)))