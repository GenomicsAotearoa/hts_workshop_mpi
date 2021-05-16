## Glossary

#### Absolute path

* A [path](#path) that refers to a particular location in a file system.
* Absolute paths are usually written with respect to the file system's [root directory](#root-directory), and begin with `/` (`\\` on Microsoft Windows).
* *See also: [relative path](#relative-path)*

#### Argument

* A value given to a function or program when it runs.
* The term is often used interchangeably (and inconsistently) with [parameter](#parameter).

#### Command shell

* See [shell](#shell)

#### Command-line interface

* A user interface based on typing commands, usually at a [REPL](#read-evaluate-print-loop).
* *See also: [graphical user interface](#graphical-user-interface)*

#### Comment

* A remark in a program that is intended to help human readers understand what is going on, but is ignored by the computer.
* Comments in `Python`, `R`, and the Unix shell start with a `#` character and run to the end of the line, other languages have other conventions.

#### Current working directory

* The directory that [relative paths](#relative-path) are calculated from.
* Equivalently, the place where files referenced by name only are searched for.
* Every [process](#process) has a current working directory.
* The current working directory is usually referred to using the shorthand notation `.` (pronounced "dot").

#### File system

* A set of files, directories, and I/O devices (such as keyboards and screens).
* A file system may be spread across many physical devices, or many file systems may be stored on a single physical device; the [operating system](#operating-system) manages access.

#### Filename extension

* The portion of a file's name that comes after the final "." character.
* By convention this identifies the file's type.
  * For example, `.txt` means "text file", `.png` means "Portable Network Graphics file"
* These conventions are not enforced by most operating systems and it is perfectly possible (but confusing!) to name an MP3 sound file `homepage.html`.
  * Many applications use filename extensions to identify the [MIME type](#mime-type) of the file, so misnaming files may cause those applications to fail.

#### Filter

* A program that transforms a stream of data.
* Many Unix command-line tools are written as filters - they read data from [standard input](#standard-input), process it, and write the result to [standard output](#standard-output).

#### Flag

* A terse way to specify an option or setting to a command-line program.
* Conventions for flags vary between operating systems.
  * Unix applications use a dash followed by a single letter, such as `-v`, or two dashes followed by a word, such as `--verbose`.
  * DOS applications use a slash, such as `/V`.
* Depending on the application, a flag may be followed by a single argument, as in `-o /tmp/output.txt`.

#### For loop

* A loop that is executed once for each value in some kind of set, list, or range.
* *See also: [loop](#loop), [while loop](#while-loop)*

#### Graphical user interface


* A user interface based on selecting items and actions from a graphical display, usually controlled by using a mouse.
* *See also: [command-line interface](#command-line-interface)*

#### Home directory

* The default directory associated with an account on a computer system.
* By convention, all of a user's files are stored in or below her home directory.

#### Loop

* A set of instructions to be executed multiple times. Consists of a [loop body](#loop-body) and (usually) a condition for exiting the loop.
* *See also [for loop](#for-loop) and [while loop](#while-loop)*

#### Loop body

* The set of statements or commands that are repeated inside a [for loop](#for-loop) or [while loop](#while-loop).

#### MIME type

*MIME (Multi-Purpose Internet Mail Extensions) types describe different file types for exchange on the Internet, for example images, audio, and documents.

#### Operating system

* Software that manages interactions between users, hardware, and software [processes](#process). Common examples are Linux, OS X, and Windows.

#### Parameter

* A variable named in a function's declaration that is used to hold a value passed into the call.
* The term is often used interchangeably (and inconsistently) with [argument](#argument).

#### Parent directory

* The directory that "contains" the one in question.
* Every directory in a file system except the [root directory](#root-directory) has a parent.
* A directory's parent is usually referred to using the shorthand notation `..` (pronounced "dot dot").

#### Path

* A description that specifies the location of a file or directory within a [file system](#file-system).
* *See also: [absolute path](#absolute-path), [relative path](#relative-path)*

#### Pipe

* A connection from the output of one program to the input of another.
* When two or more programs are connected in this way, they are called a "pipeline".

#### Process

* A running instance of a program, containing code, variable values, open files and network connections, and so on.
* Processes are the "actors" that the [operating system](#operating-system) manages.
* The [operating system](#operating-system) typically runs many process at once, allowing each to run for a few milliseconds at a time to give the impression that they are executing simultaneously.

#### Prompt

* A character or characters display by a [REPL](#read-evaluate-print-loop) to show that it is waiting for its next command.

#### Quoting

* Using quotation marks of various kinds to prevent the shell from interpreting special characters.
  * For example, to pass the string `*.txt` to a program, it is usually necessary to write it as `'*.txt'` so that the shell will not try to expand the `*` wildcard.

#### Read-evaluate-print loop

* (REPL): A [command-line interface](#command-line-interface) that reads a command from the user, executes it, prints the result, and waits for another command.

#### Redirect

* To send a command's output to a file rather than to the screen or another command, or equivalently to read a command's input from a file.

#### Regular expression

* A pattern that specifies a set of character strings.
* REs are most often used to find sequences of characters in strings.

#### Relative path

* A [path](#path) that specifies the location of a file or directory with respect to the [current working directory](#current-working-directory).
* Any path that does not begin with a separator character (`/` or `\\`) is a relative path.
* *See also: [absolute path](#absolute-path)*

#### Root directory

* The top-most directory in a [file system](#file-system).
* Its name is `/` on Unix (including Linux and Mac OS X) and `\\` on Microsoft Windows.

#### Shell

* A [command-line interface](#command-line interface) such as Bash (the Bourne-Again Shell) or the Microsoft Windows DOS shell that allows a user to interact with the [operating system](#operating-system).

#### Shell script

* A set of [shell](#shell) commands stored in a file for re-use.
* A shell script is a program executed by the shell; the name "script" is used for historical reasons.

#### Standard input

* A process's default input stream.
  * In interactive command-line applications, it is typically connected to the keyboard.
  * In a [pipe](#pipe), it receives data from the [standard output](#standard-output) of the preceding process.

#### Standard output

* A process's default output stream.
  * In interactive command-line applications, data sent to standard output is displayed on the screen.
  * In a [pipe](#pipe), it is passed to the [standard input](#standard-input) of the next process.

#### Sub-directory

* A directory contained within another directory.

#### Tab completion

* A feature provided by many interactive systems in which pressing the Tab key triggers automatic completion of the current word or command.

#### Variable

* A name in a program that is associated with a value or a collection of values.

#### While loop

* A loop that keeps executing as long as some condition is true.
* *See also: [loop](#loop), [for loop](#for-loop)*

#### Wildcard

* A character used in pattern matching.
  * In the Unix shell, the wildcard `*` matches zero or more characters, so that `*.txt` matches all files whose names end in `.txt`.

---

## External references

### Opening a terminal

* [Using a UNIX/Linux emulator (Cygwin) or Secure Shell (SSH) client (Putty)](http://faculty.smu.edu/reynolds/unixtut/windows.html)
* [Addressing the digital divide in contemporary biology: Lessons from teaching UNIX](http://www.biorxiv.org/content/early/2017/04/07/122424.full.pdf+html)
* [Unix cheat sheet](https://files.fosswire.com/2007/08/fwunixref.pdf)

### Manuals

* [GNU BASH reference](https://www.gnu.org/software/bash/manual/html_node/index.html)
* [GNU manual](http://www.gnu.org/manual/manual.html)
* [Core GNU utilities](http://www.gnu.org/software/coreutils/manual/coreutils.html)

### FASTQ files

* [Phred quality score](https://en.wikipedia.org/wiki/Phred_quality_score)
* [FASTQ format](https://en.wikipedia.org/wiki/FASTQ_format)
