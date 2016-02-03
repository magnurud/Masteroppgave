%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[]
%{}
% This is a LaTeX document class designed to format documents of type
% FFI/RAPPORT, FFI/NOTAT and FFI/REISERAPPORT, according to the
% directions given in
% ``Bestemmelser om publikasjoner for Forsvarets Forskningsinstitutt''. 
% 
% The class can handle documents written in Norwegian or English. It is
% based on the standart 'article' document class of LaTeX.
%
% The package tabularx is used in this package in order to get figure
% captions to behave correctly.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Written and maintained by Kjetil Mosesen, Aanund Storhaug, Bjørn
% Jervell Hansen and Erlend Andreas Garberg
%
% Based on ffiarticle.cls, written by Trond Varslot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Class options
% =============
% 
% The class is invoked by
% 
% \documentclass[<class options>]{ffiarticle}
% 
% where class options can be 
% one of:
%   norsk               - norwegian language used in document
%   nynorsk             - new norwegian language used in document
%   english             - english language used in document
%
% one of:
%   rapport             - FFI/RAPPORT
%   notat               - FFI/NOTAT
%   reise               - FFI/REISERAPPORT
% 
% one of:               The paper is classified as:
%   begrenset           - Restricted
%   konfidensielt       - Confidential
%   hemmelig            - Secret
%   unntattOffentlighet - Exempt from public disclosure
%   fortrolig           - in confidence
%
% and zero or more of:
%   twoside             - Double sided printing
%   venstrejustert      - Raggedright text style, as opposed to block
%                         style
%   blokkjustert        - Block (justified) text style
%
% If no options are selected explicitly, the option 'rapport' is used,
% and the report is not classified. Default language is norwegian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{ffiarticle}[2012/08/23 FFI report class]

%%%%%%%%%%%%%%%%%%%%%%%%%
% Parametre til klassen %
%%%%%%%%%%%%%%%%%%%%%%%%%

%Standard dokumenttype er en norsk, ugradert rapport med rett
%venstremarg som godkjennes av avdelingssjef. %AFS
\def\@pubtype{rapport}
\def\@rtype{ugradert}
\def\@spraak{norsk}
\def\@tekstjustering{blokk}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definisjon av hvilke parametre dokumentklassen kan ta. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\DeclareOption{norsk}{\def\@spraak{norsk}}
\DeclareOption{nynorsk}{\def\@spraak{nynorsk}}
\DeclareOption{english}{\def\@spraak{british}}

\DeclareOption{rapport}{\def\@pubtype{rapport}}
\DeclareOption{report}{\def\@pubtype{rapport}}
\DeclareOption{notat}{\def\@pubtype{notat}}
\DeclareOption{reise}{\def\@pubtype{reiserapport}}
\DeclareOption{reiserapport}{\def\@pubtype{reiserapport}}

\DeclareOption{begrenset}{\def\@rtype{begrenset}}
\DeclareOption{konfidensielt}{\def\@rtype{konfidensielt}}
\DeclareOption{hemmelig}{\def\@rtype{hemmelig}}
\DeclareOption{fortrolig}{\def\@rtype{fortrolig}}
\DeclareOption{unntattOffentlighet}{\def\@rtype{unntattOffentlighet}}


%Case list for option FFI@downgradyears
\newcommand*{\FFI@ufan}{ufan}
\newcommand*{\FFI@twoyears}{2years}
\newcommand*{\FFI@fiveyears}{5years}
\newcommand*{\FFI@thirtyyears}{30years}

%Set time of downgrad
%TODO: To set only constants once, find a command that process command
%given in declaration
\newcommand*{\FFI@downgradyears}{30years}%{\FFI@thirtyyears}%Default value
\DeclareOption{ufan}{\renewcommand*{\FFI@downgradyears}{ufan}}%{\FFI@ufan}}
\DeclareOption{2år}{\renewcommand*{\FFI@downgradyears}{2years}}%{\FFI@twoyears}}
\DeclareOption{5år}{\renewcommand*{\FFI@downgradyears}{5years}}%{\FFI@fiveyears}}
\DeclareOption{30år}{\renewcommand*{\FFI@downgradyears}{30years}}%{\FFI@thirtyyears}}

% Text style, ragged right margin or justified
\DeclareOption{venstrejustert}{\def\@tekstjustering{venstre}}
\DeclareOption{venstre}{\def\@tekstjustering{venstre}}
\DeclareOption{blokkjustert}{\def\@tekstjustering{blokk}}
\DeclareOption{blokk}{\def\@tekstjustering{blokk}}


%La article ta seg av alle parametre som ikke er definert her.
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions
%Baserer seg paa standard article.cls
\LoadClass[11pt,a4paper]{article}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definisjon av konstanter %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\def\norsk{norsk}
\def\nynorsk{nynorsk}
\def\british{british}
\def\rapport{rapport}
\def\reiserapport{reiserapport}
\def\publication{publication}
\def\notat{notat}
\def\unntattOffentlighet{unntattOffentlighet}
\def\fortrolig{fortrolig}
\def\privat{privat}
\def\ugradert{ugradert}
\def\begrenset{begrenset}
\def\konfidensielt{konfidensielt}
\def\hemmelig{hemmelig}
\def\venstre{venstre}
\def\blokk{blokk}
\def\avdelingssjef{avdelingssjef}
\def\forskningssjef{forskningssjef}


 
%%%%%%%%%%%%%%%%%%%%%%%%
% Innlesing av pakker  %
%%%%%%%%%%%%%%%%%%%%%%%%

%Babel is used to handle different languages
\RequirePackage[\@spraak]{babel}

\RequirePackage{ifthen}
\RequirePackage{ifpdf}
\RequirePackage{calc}
\RequirePackage{fp}%floating point arithmetic

%Topp og bunntekst bruker denne pakken.
\RequirePackage{fancyhdr}


\RequirePackage{lastpage}%for referencing to last page

%Font packages
\RequirePackage{times}

\RequirePackage[titles]{tocloft}% FOr customising TOC


% For customising the space around section commands
%\RequirePackage{titlesec}

%Enable coloring of text, driver spec necessary?
\RequirePackage{color}
%Input encoding
\RequirePackage[latin1]{inputenc}

\ifx\@tekstjustering\venstre
	%\raggedright with hyphenations
	\RequirePackage{ragged2e}
	\setlength{\RaggedRightRightskip}{0pt plus 3em} %default: {0pt plus 2em}
\fi
%Microtypographic extension to pdfTeX
\ifpdf
	\RequirePackage[babel=true]{microtype}
\else
	\relax
\fi

%%%%%%%%%%%%%%%%
%
% Egne makroer 
%
%%%%%%%%%%%%%%%%
\newcommand{\versiontext}[1]{%DEPRICATED
  \def\@versiontext{{#1}}}
\newcommand{\reportnumber}[1]{
  \def\@reportno{#1}}
\newcommand{\reference}[3]{%Only used in \projectnumber
  \def\@dept{#1}
  \def\@project{#2}
  \def\@archiveno{#3}}
%En enklere kommando for kun prosjektnummer
\newcommand{\projectnumber}[1]{%
 \reference{}{#1}{}%
}

\newcommand{\englishtitle}[1]{ %DEPRICATED 
  \ifx\@rtype\ugradert
    \def\@etitle{#1}
  \else
    \ifx\@rtype\unntattOffentlighet
      \def\@etitle{#1}
    \else
      \ifx\@rtype\fortrolig
        \def\@etitle{#1}
      \else
        \def\@etitle{\@engtitleprefix #1}
      \fi %fortrolig
    \fi %unntattOffentlighet
  \fi %ugradert
}

\newcommand{\offentlighetsloven}[1]{
  \def\@offloven{#1}}

\newcommand{\subtitle}[1]{ %DEPRICATED
  \def\@subtitle{#1}} 

\newcommand{\projectleader}[1]{
  \def\@projectleader{#1}}
\newcommand{\chiefscientist}[1]{
  \def\@chiefscientist{#1}}
\newcommand{\divisionmanager}[1]{ 
  \def\@divisionmanager{#1}}      
\newcommand{\directorgeneral}[1]{ 
  \def\@directorgeneral{#1}}            
\newcommand{\dateofpublishingNo}[1]{
  \def\@dateofpublishingNo{#1}}
\newcommand{\dateofpublishingEn}[1]{
  \def\@dateofpublishingEn{#1}}
\newcommand{\downgrad}[1]{%DEPRICATED. Downgrading schedule is set as
                          %an option
\newcommand{\eref}[1]{(\ref{#1})}
  \def\@downgrad{#1}}
\newcommand{\thesaurusreference}[1]{%DEPRICATED
  \def\@thesaurusref{#1}}

% Lagt til en valgfri parameter for å fange opp egen gradering for
% sammendragene 
\newcommand{\abstractdocpage}[2][]{
  \def\@classificationSummary{#1}
  \def\@abstractdocpage{#2}}
\renewcommand{\abstract}[2][]{
  \def\@classificationSummary{#1}
  \def\@abstractdocpage{#2}}
\newcommand{\sammendrag}[2][]{
  \def\@graderingSammendrag{#1}
  \def\@sammendrag{#2}}
\newcommand{\forord}[1]{
  \def\@forord{#1}}
\newcommand{\isbn}[1]{ % Lagt til av Duc V. Duong
  \def\@isbn{#1}}
\newcommand{\isbnP}[1]{
  \def\@isbn{#1}}
\newcommand{\isbnE}[1]{
  \def\@isbnE{#1}}

\newcommand\@emptypagetext{ }
\newcommand{\startnewrightpage}%
% Makro for å sette inn en blank side (brukt ved tosidig rapport)
{
  \clearpage
  \if@twoside
    \ifodd\c@page
      {}
    \else
      \thispagestyle{empty}
      \parbox[c][\textheight]{\textwidth}{\centerline\@emptypagetext}
      \newpage
    \fi
  \fi
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% TEKST-DEFINISJONER
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Disse skal være de samme på engelsk og norsk! 

\def\@begr{BEGRENSET}
\def\@konf{KONFIDENSIELT}
\def\@hemm{HEMMELIG}
\def\@fortr{FORTROLIG}
\def\@uoff{Unntatt offentlighet}
\def\@unclassifiedtitleprefix{(U)}
\def\@restrictedtitleprefix{(B)}



\ifx\@spraak\norsk
  \def\@ugrad{}
  \def\@ffi{Forsvarets forskningsinstitutt (FFI)}
  \def\@date{\@dateofpublishingNo}
  \def\@forordtekst{Forord}
  \def\@emne{Emneord}
  \def\@godkjentAv{Godkjent av}
\else
  \ifx\@spraak\nynorsk
    \def\@ugrad{}
    \def\@ffi{Forsvarets forskningsinstitutt (FFI)}
    \def\@date{\@dateofpublishingNo}
    \def\@forordtekst{Forord}
    \def\@emne{Emneord}
    \def\@godkjentAv{Godkjend av}
  \else %engelsk
    %\renewcommand{\bibname}{KLASJDH**ASDKLASJD}
    \def\@ugrad{}
    \def\@ffi{Norwegian Defence Research Establishment (FFI)}
    \def\@date{\@dateofpublishingEn}
    \def\@forordtekst{Preface}
    \def\@emne{Keywords}
    \def\@godkjentAv{Approved by}
  \fi
\fi

% Setting the appropriate document type text
\def\@pubtypetext{}
\ifx\@pubtype\notat
  \def\@pubtypetext{notat}
\else 
  \ifx\@pubtype\rapport
    \def\@pubtypetext{rapport}
  \else
    \ifx\@pubtype\reiserapport
      \def\@pubtypetext{reiserapport}
    \fi
  \fi
\fi



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SETTING AV GRADERINGSTEKST
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\def\@emptylawquote{}
\def\@offloven{IKKE SATT} %Default-verdi for å unngå en kryptisk
                          %feilmelding hvis \offentlighetsloven ikke
                          %er satt av brukeren

%\def\@lawquote{\mbox{iht. sikkerhetsloven \S 11 og \S 12}\newline
%\mbox{jf. offentleglova \S 13}}
\def\@lawquote{\parbox{4.5cm}{\mbox{iht. sikkerhetsloven \S\ 11 og \S\ 12}
    \newline \mbox{jf. offentleglova \S\ 13}}}


\ifx\@rtype\fortrolig
  \def\@lawquote{\parbox{4.1cm}{\mbox{iht. Beskyttelsesinstruksen} \newline 
  	\mbox{jf. offentleglova \@offloven}}} 
\fi
\ifx\@rtype\unntattOffentlighet 
  \def\@lawquote{\parbox{4.3cm}{iht. offentleglova \@offloven}}
\fi


%Downgrad texts
\newcommand*{\FFI@downgradyearstext}{}%Default value (30 years)
\ifx\FFI@downgradyears\FFI@twoyears
 \renewcommand*{\FFI@downgradyearstext}{\mbox{Avgraderes etter 2 år}}
\fi
\ifx\FFI@downgradyears\FFI@fiveyears
 \renewcommand*{\FFI@downgradyearstext}{\mbox{Avgraderes etter 5 år}}
\fi
\ifx\FFI@downgradyears\FFI@ufan
 \ifx\@rtype\fortrolig
   \renewcommand*{\FFI@downgradyearstext}{
     \mbox{Unntatt fra automatisk}\newline
     \mbox{nedgradering iht.}\newline
     \mbox{Beskyttelsesinstruksen \S\ 5a}}
 \else
   \renewcommand*{\FFI@downgradyearstext}{
     \mbox{Unntatt fra automatisk}\newline
     \mbox{nedgradering iht. sikkerhetsloven }\newline
     \mbox{\S\ 11}}
 \fi
\fi



%downgrad er ikke definert her

\ifx\@rtype\ugradert 
  \def\@gradtext{\@ugrad}
  \def\@lawtext{\@emptylawquote}
\fi

\ifx\@rtype\begrenset 
  \def\@gradtext{\@begr} 
  \def\@lawtext{\@lawquote}
\fi

\ifx\@rtype\konfidensielt 
  \def\@gradtext{\@konf} 
 \def\@lawtext{\@lawquote} 
\fi

\ifx\@rtype\hemmelig 
  \def\@gradtext{\@hemm} 
  \def\@lawtext{\@lawquote} 
\fi

\ifx\@rtype\fortrolig 
  \def\@gradtext{\@fortr} 
  \def\@lawtext{\@lawquote}
\fi

\ifx\@rtype\unntattOffentlighet 
  \def\@gradtext{\@uoff}
  \def\@lawtext{\@lawquote}
\fi 



% Kommandoer for å sette tittelen til ugradert (U) eller begrenset (B)

\def\@titleprefix{\@empty}
\def\@engtitleprefix{\@empty}

\newcommand{\restrictedtitle}{
  \def\@titleprefix{\@restrictedtitleprefix}
  \def\@engtitleprefix{\@engrestrictedtitleprefix}}
\newcommand{\unclassifiedtitle}{
  \def\@titleprefix{\@unclassifiedtitleprefix}
  \def\@engtitleprefix{\@unclassifiedtitleprefix}}


%%%%%%%%%%%%%%%%
%
% FAKTABOKS
%
%%%%%%%%%%%%%%%%


\newcommand{\faktaboks}[1]{
\par\hspace*{-9mm}
\setlength{\fboxrule}{0.5pt}
\setlength{\fboxsep}{6mm}
\fcolorbox[gray]{0.9}{0.875}{
\begin{minipage}{\linewidth}
\color{ffiDarkBlue}{#1}
\end{minipage}
}}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EMNEORD
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand{\key}[2]{% Only used in \emneord
  \stepcounter{refcounter}
  \addcontentsline{engkey}{refkey}{#1}%
  \addcontentsline{norkey}{refkey}{#2}%
}

\newcommand{\emneord}[1]{%
 \key{}{#1}
}


\newcommand{\l@refkey}[2]{#1\par\vspace*{0.5cm}}
\newcounter{refcounter}
\newcommand{\@norkey}{\@starttoc{norkey}}
\newcommand{\@engkey}{\@starttoc{engkey}}



%%%%%%%%%%%%%%%%%%%%%%%%
% Forside              %
%%%%%%%%%%%%%%%%%%%%%%%%


\def\@forside{
  \maketitle
}


% Distances
\newcommand{\headerToTitle}{3.8cm}
\newcommand{\titleToAuthor}{1.9cm}
\newcommand{\authorToFFI}{9.8cm}
\newcommand{\FFIToDate}{0.2cm}

% Definition of title page
\renewcommand\@maketitle{
  \vfill
  \begin{center}
    \null
    \vskip \headerToTitle
    {\color{ffiDarkBlue}\sffamily\Large\bfseries \@titleprefix\ \@title\par}%
    \vskip \titleToAuthor
    {\@author\par}%
    \vskip \authorToFFI
    {\@ffi\par}%
    \vskip \FFIToDate
    {\@date\par}%
  \end{center}
  \vfill
}

% Redefine \maketitle to get correct header and footer
\renewcommand\maketitle{%
  \def\@ufanquote{Unntatt fra automatisk nedgradering\\iht. sikkerhetsloven § 11}
  \def\@twoyearquote{Avgraderes etter 2 år}
  \def\@fiveyearquote{Avgraderes etter 5 år}
  \def\@thirtyyearquote{}
  
  \def\@emptydowngradquote{}
  \def\@downgradquote{}
  
  \ifx\@downgrad\ufan \def\@downgradquote{\@ufanquote} \fi
  \ifx\@downgrad\twoyears \def\@downgradquote{\@twoyearquote} \fi
  \ifx\@downgrad\fiveyears \def\@downgradquote{\@fiveyearquote} \fi
  \ifx\@downgrad\thirtyyears \def\@downgradquote{\@emptydowngradquote} \fi
  
  % Only use downgrad quote if document is classified
  \ifx\@rtype\ugradert \def\@downgradtext{\@emptydowngradquote} \fi
  \ifx\@rtype\begrenset \def\@downgradtext{\@downgradquote} \fi
  \ifx\@rtype\konfidensielt \def\@downgradtext{\@downgradquote} \fi
  \ifx\@rtype\hemmelig \def\@downgradtext{\@downgradquote} \fi
  \ifx\@rtype\fortrolig \def\@downgradtext{\@emptydowngradquote} \fi
  \ifx\@rtype\unntattOffentlighet \def\@downgradtext{\@emptydowngradquote} \fi

  \fancypagestyle{plain}{\ffiTitleHeaderFooter}%
  \thispagestyle{plain}%
  \@maketitle%
}

%%%%%%%%%%%%%%%%%%%%%%%%
%
% Direktørens side     
%
%%%%%%%%%%%%%%%%%%%%%%%%

% Distances
\newcommand{\headerToReport}{0.1cm}
\newcommand{\reportToProject}{0.4cm}
\newcommand{\projectToIsbnP}{0.8cm}
\newcommand{\isbnPToIsbnE}{-0.4cm}
\newcommand{\isbnEToKeywords}{1.8cm}
\newcommand{\keywordpre}{0.4cm}
\newcommand{\keywordsize}{7.1cm}
\newcommand{\linedistance}{0.5cm}

% Define empty value for ISBN
\def\@isbn{}  % Lagt til av Duc V. Duong
\def\@isbnE{}

% Define tomme verdier  for Approved by
\def\@projectleader{}
\def\@chiefscientist{}
\def\@divisionmanager{}
\def\@directorgeneral{}
\def\@testus{}

\def\@dirside{
  \newpage
  \null
  \vskip \headerToReport%
  {FFI-\@pubtypetext\ \@reportno\par}%
  \vskip \reportToProject
  {\@project\par}%
  \vskip \projectToIsbnP
  \ifx\@rtype\ugradert  % Lagt til av Duc V. Duong
    \ifx\@pubtype\rapport
      {P: ISBN\ \@isbn\par}%
	\vspace \isbnPToIsbnE
      {E: ISBN\ \@isbnE}%
    \else
      {\vskip 4em \par}%
    \fi
  \else
    {\vskip 4em \par}%
  \fi
  \vskip \isbnEToKeywords

 \ifx\@pubtype\reiserapport   
   \ifthenelse{\value{refcounter} > 0}
    {{\color{ffiDarkBlue}\sffamily\Large\textbf{\@emne}}}
    {}
  \else
    {\color{ffiDarkBlue}\sffamily\Large\textbf{\@emne}}%
  \fi

  \vskip \keywordpre
  
  \begin{minipage}[t][\keywordsize][t]{10cm}
    \@norkey
  \end{minipage}
  
  {\color{ffiDarkBlue}\sffamily\Large \textbf{\@godkjentAv}}%
  \vskip 0.3cm

  \parbox[t][5cm][t]{6.9cm}{
    \ifx\@projectleader\@empty
      \relax
    \else
      \@projectleader\par
      \vskip \linedistance
    \fi

    \ifx\@chiefscientist\@empty
      \relax
    \else
      \@chiefscientist\par
      \vskip \linedistance
    \fi

    \ifx\@divisionmanager\@empty
      \relax
    \else
      \@divisionmanager\par
      \vskip \linedistance
    \fi

    \ifx\@directorgeneral\@empty
      \relax
    \else
      \@directorgeneral\par
    \fi
  }
  \parbox[t][5cm][t]{6.9cm}{
    \ifx\@projectleader\@empty
      \relax
    \else
      \ifx\@spraak\norsk
        Prosjektleder\par
      \else 
        \ifx\@spraak\nynorsk
          Prosjektleiar\par
        \else%engelsk
          Project manager\par        
        \fi
      \fi
      \vskip \linedistance
    \fi

    \ifx\@chiefscientist\@empty
      \relax
    \else
      \ifx\@spraak\norsk
        Forskningssjef\par
      \else 
        \ifx\@spraak\nynorsk
          Forskingssjef\par
        \else%engelsk
          Director of Research\par
        \fi
      \fi
      \vskip \linedistance
    \fi

    \ifx\@divisionmanager\@empty
      \relax
    \else
      \ifx\@spraak\british
        Director\par
      \else%norsk eller nynorsk
        Avdelingssjef\par
      \fi
      \vskip \linedistance
    \fi

    \ifx\@directorgeneral\@empty
      \relax
    \else
      \ifx\@spraak\british
        Director General\par
      \else %norsk eller nynorsk
        Adm.\ direktør\par
      \fi
    \fi
  }
}%dirside


%%%%%%%%%%%%%%%%%%%%%%%%
%
% Norsk sammendrag     
%
%%%%%%%%%%%%%%%%%%%%%%%%
\def\@norskSammendrag{
 \ifx\@sammendrag\@empty
   \newpage
  \else
    \ifx\@spraak\nynorsk
      \ifx\@graderingSammendrag\@empty
	\section*{Samandrag}
      \else
        \section*{(\@graderingSammendrag ) Samandrag}
      \fi
    \else%bokmaal
      \ifx\@graderingSammendrag\@empty
	\section*{Sammendrag}
      \else
        \section*{(\@graderingSammendrag ) Sammendrag}
      \fi
    \fi
    \@sammendrag
    \newpage
  \fi
}

%%%%%%%%%%%%%%%%%%%%%%%%
%
% Engelsk sammendrag   
%
%%%%%%%%%%%%%%%%%%%%%%%%
\def\@engelskSammendrag{
  \ifx\@abstractdocpage\@empty
    \newpage
  \else
    \ifx\@classificationSummary\@empty
	\section*{English summary}
    \else
	\section*{(\@classificationSummary ) English summary}
    \fi
    \@abstractdocpage  
    \newpage
  \fi
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Forord/preface (valgfritt) 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\def\@forord{} % For å klare å ikke ha med forord dersom det ikke
                                % er med i kildefilene
\def\@preface{
  \ifx\@forord\@empty
    \relax
  \else
    \newpage
    \section*{\@forordtekst}
    \addcontentsline{toc}{section}{\protect\numberline{\null}{\@forordtekst}}
    \@forord
  \fi
}

%%%%%%%%%%%%%%%%%%%%%%%%
%
% Innholdsfortegnelse  
%
%%%%%%%%%%%%%%%%%%%%%%%%


%Parameters specifying the visual appearance of all entries
\setcounter{tocdepth}{3}%insert 3 levels, section -> subsubsection
\renewcommand\@dotsep{1000}%separation between dots - big number make dots dissappear!



%The tableofcontents!

\renewcommand\tableofcontents{
  \cleardoublepage%Remove naming of table of contents!
  \section*{\contentsname %
    \@mkboth{%
      \MakeUppercase\contentsname}{\MakeUppercase\contentsname}}%
  \vspace{1.2cm}
  \newlength{\oldparskip}
  \setlength{\oldparskip}{\the\parskip} % A hack to make toc independent
  \setlength{\parskip}{0mm}             % of the global parskip


  \renewcommand\cftdotsep{1000}

  % Avstand før hver entry:
  \setlength{\cftbeforesecskip}{4mm}
  \setlength{\cftbeforesubsecskip}{1mm}
  \setlength{\cftbeforesubsubsecskip}{1mm}
  \setlength{\cftbeforeparaskip}{1mm}

  % Avstand mellom kapittelnummer og kapitteltittel:
  \setlength{\cftsecindent}{0mm}
  \setlength{\cftsubsecindent}{0mm}
  \setlength{\cftsubsubsecindent}{0mm}
  \setlength{\cftparaindent}{0mm}

  % Boks holdt av til kapittelnummer
  \setlength{\cftsecnumwidth}{27mm}
  \setlength{\cftsubsecnumwidth}{27mm}
  \setlength{\cftsubsubsecnumwidth}{27mm}
  \setlength{\cftparanumwidth}{27mm}

  % Fonten til entriene:
  %\renewcommand\cftsecfont{\color{ffiDarkBlue}\large\sffamily\bfseries}
  \renewcommand\cftsecfont{\color{ffiDarkBlue}\fontsize{12}{15}\selectfont\sffamily\bfseries}
  %\renewcommand\cftsecpagefont{\color{ffiDarkBlue}\large\sffamily\bfseries}
  \renewcommand\cftsecpagefont{\color{ffiDarkBlue}\fontsize{12}{15}\selectfont\sffamily\bfseries}
  \renewcommand\cftsubsecfont{\fontsize{10}{15}\normalfont\sffamily}
  \renewcommand\cftsubsecpagefont{\fontsize{10}{15}\normalfont\sffamily}
  \renewcommand\cftsubsubsecfont{\fontsize{10}{15}\normalfont\sffamily}
  \renewcommand\cftsubsubsecpagefont{\fontsize{10}{15}\normalfont\sffamily}
  \renewcommand\cftparafont{\fontsize{10}{15}\normalfont\sffamily}
  \renewcommand\cftparapagefont{\fontsize{10}{15}\normalfont\sffamily}

  \@starttoc{toc}%

  % Setter den globale parskip tilbake til det den var før toc ble laget:
  \setlength{\parskip}{\the\oldparskip}
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Document Layout                  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=============================================================
%Page Layout

%Declare MS Word Parameters! 
%"Utskriftsformat->Marger"
\newcommand{\wordTop}{2.9cm}%Topp
\newcommand{\wordBottom}{3.4cm}%Bunn
\newcommand{\wordInner}{3.25cm}%Innvendig
\newcommand{\wordOuter}{2.5cm}%Utvendig
%"Utskriftsformat->Oppsett"
\newcommand{\wordFooter}{0.5cm}%Bunntekst
%"Spalter"
\newcommand{\wordColumnSep}{1cm}%Avstand
\newcommand{\wordColumnWidth}{7.25cm}%Bredde

%put origo for margins in upper left corner of page (default is 
%(1in,1in) from upper left corner
%\setlength{\hoffset}{-1in}
%\setlength{\voffset}{-1in}%Har en dårlig effekt på en del pakker (pdfpages)

%Vertical lengths
\setlength{\topmargin}{1cm - 1in}%upper margin to top of header

\ifx\@rtype\hemmelig
   \setlength{\headheight}{15mm}%height of head
\else
   \setlength{\headheight}{10mm}%height of head (Word->guessed;-) 
\fi


\setlength{\headsep}{\wordTop - (\headheight + \topmargin) - 1in}%distance from bottom of
                                %headline to top of body
\setlength{\topskip}{1ex}%distance from top of body to baseline of
                         %first line of text
\setlength{\textheight}{\paperheight - (\wordTop + \wordBottom)}%height of body


\setlength{\footskip}{1.7cm}%distance from bottom of
                            %body to bottom of foot (Word->measured)

%\setlength{\footskip}{10pt}%distance from bottom of
%                            %body to bottom of foot. Set to prevent
%                            %latex from complaining

%Horisontal lengths
\setlength{\textwidth}{\paperwidth - (\wordInner + \wordOuter)}%width of body
\setlength{\columnsep}{\wordColumnSep}%for twocolumns, separations of columns
\setlength{\oddsidemargin}{\wordInner - 1in}%left margin for odd pages
\setlength{\evensidemargin}{\wordOuter - 1in}%left margin for even pages

%=============================================================
%Text Layout
%Change the line spacing for whole document NB will influence ALL text!
\renewcommand{\baselinestretch}{1.2}

%Lengder for avsnitt
\setlength{\parindent}{0pt}
\setlength{\parskip}{3mm} % OBS: Husk at endring på parskip også fører
                          % til endring på avstandene rundt
                          % overskrifter (\@startsection)! (Se LaTeX
                          % companion s. 25)

\ifx\@tekstjustering\venstre
  \RaggedRight
\fi


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Section Commands
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\setcounter{secnumdepth}{4}%numbering down to level 4 (paragraph for
                           %article class)
%\@startsection{sec-name}{level}{indent}{pre-skip}{post-skip}{formatting}
%pre-skip: negative value toggle no! indentation of following paragraph
%post-skip: negative value => no line break

% Keep in mind that for both the pre and post skips, the actual skip
% around a section header is the skip defined in the \@startsection
% command PLUS the defined \parskip!! 

\renewcommand{\section}{%
\setcounter{figure}{0}%
\setcounter{table}{0}%
\setcounter{equation}{0}%
%\@startsection{section}{1}{0pt}{5mm}{0.01mm}{\normalfont\sffamily\bfseries\Large\color{ffiDarkBlue}}%
\@startsection{section}{1}{0pt}{5mm}{0.01mm}{\fontsize{14.2}{17}\selectfont\sffamily\bfseries\color{ffiDarkBlue}}%
}
 

\renewcommand{\subsection}{%
\@startsection{subsection}{2}{0pt}{3mm}{1sp}{\normalfont\sffamily\bfseries}%
} 

\renewcommand{\subsubsection}{%
\@startsection{subsubsection}{3}{0pt}{2mm}{1sp}{\normalfont\sffamily}%
} 

\renewcommand{\paragraph}{%
\@startsection{paragraph}{4}{0pt}{2mm}{1sp}{\normalfont\sffamily}%
}


%------ test for å fjerne appendix i venstrekolonna i toc
%	\newcommand{\nocontentsline}[3]{}
%	\newcommand{\tocless}[2]
%		 {\bgroup\let\addcontentsline=\nocontentsline#1{#2}\egroup}
%test:
 %\addtocontents{toc}{\hspace{27mm}}
 %\addcontentsline{toc}{section}{}
 		
%\color{ffiDarkBlue}

		
%%%
% Code to get letters for numbers in appendix section headings and
% make the top level appendix section prefixed with 'Appendix'
% while keeping the subsections, figures, and tables prefixless
%%% %\renewcommand\thesection{Appendiks \@Alph\c@section} 	
\renewcommand{\appendix}{%

	\ifx\@spraak\british
	  \renewcommand{\appendixname}{Appendix}
	\else %norsk eller nynorsk
	  \renewcommand{\appendixname}{Vedlegg}
	\fi

  %\setlength{\cftsecindent}{5cm}

	\renewcommand\thesection{\appendixname\ \@Alph\c@section}
	\renewcommand\thesubsection{\@Alph\c@section.\arabic{subsection}}
	
	\renewcommand\thefigure{\@Alph\c@section.\arabic{figure}}
	\renewcommand\thetable{\@Alph\c@section.\arabic{table}}
	\renewcommand\theequation{\@Alph\c@section.\arabic{equation}}
	
	\setcounter{section}{0} % To start at 'A'

%%\cftsetindents{section}{10cm}{10cm}

%	\makeatletter
%	\let\stdl@section\l@section
%	\renewcommand*{\l@section}[2]{
%		\l@section{\textcolor{green}{#1}}{\textcolor{blue}{#2}}
%	}
%	%	\makeatother


%%fråå ieeetrans.cls
%\def\tableofcontents{\section*{Contents}\@starttoc{toc}}
%\def\l@section#1#2{\addpenalty{\@secpenalty} \addvspace{1.0em plus 1pt}
%    \@tempdima 1.8em \begingroup \parindent \z@ \rightskip \@pnumwidth
%    \parfillskip-\@pnumwidth \bf\leavevmode #1\hfil\hbox to\@pnumwidth{\hss #2}
%    \par \endgroup}
%\def\l@subsection{\@dottedtocline{2}{1.5em}{2.3em}}
%\def\l@subsubsection{\@dottedtocline{3}{3.8em}{3.2em}}
%\def\listoffigures{\section*{List of Figures}\@starttoc{lof}}
%\def\l@figure{\@dottedtocline{1}{1em}{1.8em}}
%\def\listoftables{\section*{List of Tables}\@starttoc{lot}}
%\let\l@table\l@figure

}

	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Header and Footer declaration
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\def\@reportno{xxx/yyy}

\newcommand{\ffiHeaderFooter}{%define the header and footer
\fancyhf{}%clear all header and footer fields

\ifx\@rtype\fortrolig 
  \relax
  \else
  \ifx\@rtype\unntattOffentlighet
    \relax
    \else
    \ifx\@rtype\hemmelig
      \fancyhead[RO,LE]{\textsf{\large \bfseries \@gradtext}\\
        \footnotesize\textsf{Side {\thepage} av \pageref{LastPage}}\\}
      \fancyfoot[RO,LE]{\textsf{\large \bfseries \@gradtext}}
    \else
    \fancyhf[RO,LE]{\textsf{\large \bfseries \@gradtext}}%print
  %classification on outer page top/bottom
    \fi
  \fi
\fi

\fancyhfoffset[RO,LE]{0cm}%adjust the header and footer outer
                          %placement relative to \textwidth
\fancyfoot[C]{\textsf{\bfseries \thepage}}%print pagenumber in center of footer
\fancyfoot[RE,LO]{\textsf{\footnotesize FFI-\@pubtypetext\ 
%\fancyfoot[RE,LO]{\textsf{\fontsize{8}{12} FFI-\@pubtypetext\ 
    \@reportno}}%TODO: english text?

\renewcommand{\headrulewidth}{0pt}
\renewcommand{\footrulewidth}{0pt}
}%end ffiHeaderFooter


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Header and Footer declaration for titlepage
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newcommand{\ffiTitleHeaderFooter}{%define the header and footer
\fancyhf{}%clear all header and footer fields


%Finding color according to the restriction on the document:
\def\@titlepagecolor{ffiClassified}
\ifx\@rtype\fortrolig \def\@titlepagecolor{ffiDarkBlue} \fi
\ifx\@rtype\unntattOffentlighet \def\@titlepagecolor{ffiDarkBlue} \fi

\fancyhead[LO,RE]{
{\fontsize{8}{12}\sffamily FFI-\@pubtypetext} {\fontsize{8}{12}\sffamily \@reportno\\} % A line break is added to align the left header
                    % with the right header vertically
}

\fancyhead[RO,LE]{
  %\begin{tabular}{p{3.9cm}} % To make the text flushed left
  \begin{tabular}{p{41mm}} % To make the text flushed left
    \textsf{\large\bfseries\color{\@titlepagecolor}\@gradtext}\\
    \textsf{\footnotesize\color{\@titlepagecolor}\@lawtext}\\
%If downgrad 30 years, then no downgrad text!
\ifx\FFI@downgradyears\FFI@thirtyyears
\relax
\else
    \textsf{\footnotesize\color{\@titlepagecolor}\FFI@downgradyearstext}\\
\fi
%HEMMELIG, m.m? also adds page number and total pages on each
%page in the "classification header"
\ifx\@rtype\hemmelig
    \textsf{\footnotesize\color{\@titlepagecolor} Side {\thepage} av
      \pageref{LastPage}}
\fi
  \end{tabular}
  }
 

\ifx\@rtype\unntattOffentlighet 
  \relax
  \else
  \ifx\@rtype\fortrolig
    \relax
    \else
       \fancyfoot[RO,LE]{
\begin{tabular}{p{3.5cm}}%To left-adjust equal to header
\textsf{\large\bfseries\color{\@titlepagecolor}
             \@gradtext} %print classification on outer page
                          %top/bottom
\end{tabular}
}
  \fi
\fi



% Hack to get left justified classification on right side of header
%\fancyhfoffset[LO,RE]{-9.5cm}%adjust the header and footer outer
                          %placement relative to \textwidth

% Hack to set the left hand header piece
%\fancyhfoffset[RO,LE]{3cm}

\setlength{\headheight}{2.5cm}
\setlength{\headsep}{\wordTop - \headheight}
%\addtolength{\footskip}{-\headheight} % For å forhindre at footeren
%havner for langt nede (BHn 28.10.2008): Fjernet fordi den genererte
%feilmelding. Ser ikke ut til at footeren blir feil nå...

\renewcommand{\headrulewidth}{0pt}
\renewcommand{\footrulewidth}{0pt}
}%end ffiTitleHeaderFooter

\fancypagestyle{plain}{\ffiHeaderFooter}%redefine plain style since it
                                %is used by \maketitle

\pagestyle{plain}%toggle plain style for whole document



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Counters redefinition
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\renewcommand{\thefigure}{\thesection.\arabic{figure}}
\renewcommand{\thetable}{\thesection.\arabic{table}}
\renewcommand{\theequation}{\thesection.\arabic{equation}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Caption - figure and table text
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\RequirePackage{caption}
%\DeclareCaptionLabelSeparator{ffistyle}{\hspace{4mm}}
\ifx\@tekstjustering\venstre
  \captionsetup{justification=raggedright,font=it,format=hang,singlelinecheck=false,labelsep=quad}
\else
  \captionsetup{font=it,format=hang,singlelinecheck=false,labelsep=quad}
  % labelsep=ffistyle
\fi


%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%
% Definition of FFI colors
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%macro for integer representation of RGB colors
\newcommand{\definecolorRGB}[4]{\FPupn\RGBr{255 #2 /}\FPupn\RGBg{255
    #3 /}\FPupn\RGBb{255 #4
    /}\definecolor{#1}{rgb}{\RGBr,\RGBg,\RGBb}}

\definecolorRGB{ffiClassified}{255}{0}{0}
\definecolorRGB{ffiYellow}{243}{224}{88}
\definecolorRGB{ffiLightBlue}{122}{162}{201}
\definecolorRGB{ffiDarkBlue}{56}{89}{119}
\definecolorRGB{ffiLightGreen}{198}{195}{73}
\definecolorRGB{ffiDarkGreen}{102}{132}{125}
\definecolorRGB{ffiOrange}{194}{147}{87}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SINGLE SPACING IN LISTS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Based on package by Hydro-Quebeck

% Original enviroments
\let\orig@Itemize =\itemize
\let\orig@Enumerate =\enumerate
\let\orig@Description =\description

% Zero vertical spacing params
\def\Nospacing{\itemsep=0pt\parskip=0pt\parsep=0pt}

% Redefine environment
\renewenvironment{itemize}{\orig@Itemize\Nospacing}{\endlist}
\renewenvironment{enumerate}{\orig@Enumerate\Nospacing}{\endlist}
\renewenvironment{description}{\orig@Description\Nospacing}{\endlist}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% TITLEPAGE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newcommand{\ffititlepage}
{\@forside
\@dirside
\newpage
\ifx\@spraak\british
 \@engelskSammendrag
 \@norskSammendrag
\else%norsk eller nynorsk
 \@norskSammendrag
 \@engelskSammendrag
\fi

%\ifx\@spraak\norsk
% \@norskSammendrag
% \@engelskSammendrag
%\else
% \@engelskSammendrag
% \@norskSammendrag
%\fi
\tableofcontents
\@preface
\cleardoublepage % To make sure section 1 starts at a odd-numbered page
\vspace*{-1cm}% A hack to make sure that the first section title end
              % up in approxiamtely the same place as the first section
              % title using Word
}

%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DEFAULT VALUES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%Default values for some variables which if skipped (forgotten by author?) will cause LaTeX to fail
\sammendrag{Du må fylle inn sammendrag. I de tilfeller der det ikke
  ønskes norsk sammedrag, må sammendraget defineres som tomt:
  \texttt{$\backslash$sammendrag\{\}}.}
\abstractdocpage{Normally, an english summary is needed. An English
  summary is avoided by defining 
  $\backslash$abstract as empty: \texttt{$\backslash$abstract\{\}}.} 
\dateofpublishingNo{Fyll inn publiseringsdato}
\dateofpublishingEn{Add date of publishing}
\reportnumber{Rapportnummer?}
\projectnumber{Her skal prosjektnummer stå}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mest sannsynlig bare dummykommandoer som defineres for å få          %
% dokumenter skrevet etter den gamle malen til å kompilere også i den  %
% nye malen.                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand{\numberofpublications}[1]{} %Tidligere brukt i
                                %fordelingslista for dokumenter
                                %gradert KONFIDENSIELT eller høyere
\newcommand{\totalnumber}[1]{}
\newcommand{\restricteddocpage}{}
\newcommand{\unclassifieddocpage}{}
\newcommand{\internaldist}[2]{} % Tidligere brukt i distribusjonslista
\newcommand{\externaldist}[2]{} % Tidligere brukt i distribusjonslista
