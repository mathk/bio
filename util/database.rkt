#lang typed/racket

(provide string->dna)

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
