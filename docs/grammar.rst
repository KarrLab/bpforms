.. _grammar:

`BpForms` grammar
------------------

The `BpForms` grammar unambiguously represents the primary structure of biopolymer forms that contain canonical and non-canonical monomeric forms using (a) a syntax similar to IUPAC/IUBMB and (b) extended alphabets for DNA, RNA, and proteins to describe monomeric forms.

`BpForms.org <https://www.bpforms.org>`_ contains a detailed overview of the grammar and examples.


Grammar
^^^^^^^

The following is the definition of the `BpForms` grammar. The grammar is defined in `Lark syntax <https://lark-parser.readthedocs.io/en/latest/grammar/>`_ which is based on `EBNF syntax <https://en.wikipedia.org/wiki/Extended_Backus%E2%80%93Naur_form>`_.

.. literalinclude:: ../bpforms/grammar.lark
    :language: text
