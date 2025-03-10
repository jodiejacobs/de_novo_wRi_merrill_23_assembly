��    &      L  5   |      P  8   Q  B   �  A   �  6     H   F  I   �  F   �  9      7   Z  6   �  M   �  L     O   d  H   �  {   �     y  �   �  e   `  :   �    	  �  
  �  �     �     �  <   �  1   �  &   "     I  "   X  9   {  I   �  �   �     �     �     �     �     �  �  �  <   �  C   )  A   m  ;   �  D   �  D   0  N   u  9   �  7   �  5   6  D   l  P   �  L     D   O  w   �       �     i   �  =   Y  1  �    �  �  �     �     �  �   �  1   W  *   �     �  "   �  8   �  E     �   e     �          !     5     E              !          #                                                      %                                        
                  $      	   &                 "                  -E                        (ignored for compatibility)
   -V, --version               output version information and exit
   -V, --version             display version information and exit
   -c, --context=CONTEXT     specify context for MSGID
   -d, --domain=TEXTDOMAIN   retrieve translated message from TEXTDOMAIN
   -d, --domain=TEXTDOMAIN   retrieve translated messages from TEXTDOMAIN
   -e                        enable expansion of some escape sequences
   -h, --help                  display this help and exit
   -h, --help                display this help and exit
   -n                        suppress trailing newline
   -v, --variables             output the variables occurring in SHELL-FORMAT
   COUNT                     choose singular/plural form based on this value
   MSGID MSGID-PLURAL        translate MSGID (singular) / MSGID-PLURAL (plural)
   [TEXTDOMAIN]              retrieve translated message from TEXTDOMAIN
   [TEXTDOMAIN] MSGID        retrieve translated message corresponding
                            to MSGID from TEXTDOMAIN
 Bruno Haible Copyright (C) %s Free Software Foundation, Inc.
License GPLv3+: GNU GPL version 3 or later <%s>
This is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law.
 Display native language translation of a textual message whose grammatical
form depends on a number.
 Display native language translation of a textual message.
 If the TEXTDOMAIN parameter is not given, the domain is determined from the
environment variable TEXTDOMAIN.  If the message catalog is not found in the
regular directory, another location can be specified with the environment
variable TEXTDOMAINDIR.
Standard search directory: %s
 If the TEXTDOMAIN parameter is not given, the domain is determined from the
environment variable TEXTDOMAIN.  If the message catalog is not found in the
regular directory, another location can be specified with the environment
variable TEXTDOMAINDIR.
When used with the -s option the program behaves like the 'echo' command.
But it does not simply copy its arguments to stdout.  Instead those messages
found in the selected catalog are translated.
Standard search directory: %s
 In normal operation mode, standard input is copied to standard output,
with references to environment variables of the form $VARIABLE or ${VARIABLE}
being replaced with the corresponding values.  If a SHELL-FORMAT is given,
only those environment variables that are referenced in SHELL-FORMAT are
substituted; otherwise all environment variables references occurring in
standard input are substituted.
 Informative output:
 Operation mode:
 Report bugs in the bug tracker at <%s>
or by email to <%s>.
 Substitutes the values of environment variables.
 Try '%s --help' for more information.
 Ulrich Drepper Usage: %s [OPTION] [SHELL-FORMAT]
 Usage: %s [OPTION] [TEXTDOMAIN] MSGID MSGID-PLURAL COUNT
 Usage: %s [OPTION] [[TEXTDOMAIN] MSGID]
or:    %s [OPTION] -s [MSGID]...
 When --variables is used, standard input is ignored, and the output consists
of the environment variables that are referenced in SHELL-FORMAT, one per line.
 Written by %s.
 error while reading "%s" missing arguments standard input too many arguments Project-Id-Version: gettext-runtime 0.20.2
Report-Msgid-Bugs-To: bug-gettext@gnu.org
PO-Revision-Date: 2020-04-17 07:38+0100
Last-Translator: Pedro Albuquerque <pmra@protonmail.com>
Language-Team: Portuguese <translation-team-pt@lists.sourceforge.net>
Language: pt
MIME-Version: 1.0
Content-Type: text/plain; charset=UTF-8
Content-Transfer-Encoding: 8bit
X-Bugs: Report translation errors to the Language-Team address.
Plural-Forms: nplurals=2; plural=(n != 1);
X-Generator: Geany / PoHelper 1.36
   -E                        (ignorado para compatibilidade)
   -V, --version               mostra informação da versão e sai
   -V, --version             mostra informação da versão e sai
   -c, --context=CONTEXTO    especifica contexto para MSGID
   -d, --domain=DOMTEXTO     obtém a mensagem traduzida de DOMTEXTO
   -d, --domain=DOMTEXTO..   obtém mensagens traduzidas de DOMTEXTO
   -e                        activa expansão de algumas sequências de escape
   -h, --help                  mostra esta mensagem e sai
   -h, --help                mostra esta mensagem e sai
   -n                        suprime nova linha final
   -v, --variables             mostra as variáveis em FORMATO-SHELL
   TOTAL                     escolhe a forma singular/plural baseado neste valor
   MSGID MSGID-PLURAL        traduz MSGID (singular) / MSGID-PLURAL (plural)
   [DOMTEXTO]                obtém a mensagem traduzida de DOMTEXTO
   [DOMTEXTO] MSGID          obtém a mensagem traduzida correspondente
                            a MSGID de DOMTEXTO
 Bruno Haible Copyright (C) %s Free Software Foundation, Inc.
Licença GPLv3+: GNU GPL versão 3 ou posterior <%s>
Este software é grátis: pode alterá-lo e redistribuí-lo.
Não há QUALQUER GARANTIA, até ao limite da Lei.
 Mostrar tradução de idioma nativo de uma mensagem textual cuja forma
gramatical depende de um número.
 Mostrar tradução de idioma nativo de uma mensagem textual.
 Se o parâmetro DIALECTO não for fornecido, o domínio é determinado a partir
da variável de ambiente TEXTDOMAIN. Se o catálogo de mensagens não for
encontrado na pasta habitual, pode ser especificada outra localização através
da variável de ambiente TEXTDOMAINDIR.
Pasta padrão de pesquisa: %s
 Se o parâmetro DIALECTO não for fornecido, o domínio é determinado a
partir da variável de ambiente TEXTDOMAIN. Se o catálogo de mensagens
não for encontrado na pasta habitual, pode ser especificada outra
localização através da variável de ambiente TEXTDOMAINDIR.
Quando usado com a opção -s, o programa comporta-se como o comando "echo".
Mas ele não copia simplesmente os seus argumentos para a saída padrão.
Ao invés, as mensagens encontradas no catálogo seleccionado são traduzidas.
Pasta padrão de pesquisa: %s
 Em modo de operação normal, a entrada padrão é copiada para a saída padrão,
com referências a variáveis de ambiente da forma $VARIÁVEL ou ${VARIÁVEL}
sendo substituídas com os valores correspondentes. Se um FORMATO-SHELL for
dado, apenas essas variáveis de ambiente referenciadas em FORMATO-SHELL são
substituídas; caso contrário, todas as variáveis de ambiente referenciadas
ocorrentes na entrada padrão são substituídas.
 Saída informativa:
 Modo de operação:
 Reporte erros no sistema de rastreio em <%s>
ou por email para <%s>.
Reporte erros de tradução em <translation-team-pt@lists.sourceforge.net>
 Substitui os valores das variáveis de ambiente.
 Tente "%s --help" para mais informação.
 Ulrich Drepper Uso: %s [OPÇÃO] [FORMATO-SHELL]
 Uso: %s [OPÇÃO] [DIALECTO] MSGID MSGID-PLURAL NÚMERO
 Uso: %s [OPÇÃO] [[DIALECTO] MSGID]
ou:  %s [OPÇÃO] -s [MSGID]...
 Quando --variables é usado, a entrada padrão é ignorada e a saída consiste
nas variáveis de ambiente referenciadas em FORMATO-SHELL, uma por linha.
 Escrito por %s.
 erro ao ler "%s" argumentos em falta entrada padrão demasiados argumentos 