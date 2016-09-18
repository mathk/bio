#lang typed/racket

(provide string->dna
         strings->dnas)

(require racket/format 
         "type.rkt")

(define (string->dna [text : String]) : Dna
  (list->vector
   (map (lambda [a] : DnaAtom
          (case a
            [(#\A) 'A]
            [(#\T) 'T]
            [(#\C) 'C]
            [(#\G) 'G]
            [else (error "Not a dna alphabet " a)]))
        (remove* '(#\newline) (string->list text)))))

(define strings->dnas
  (case-lambda
    [([textes : (Listof String)])
     (map string->dna textes)]
    [(textes : String *)
     (map string->dna textes)]))
