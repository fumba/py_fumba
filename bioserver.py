#!/usr/local/bin/python

from calculator import *
import cgi
import cgitb
cgitb.enable()


def main():
    # extract the form from the user

    form = cgi.FieldStorage()

    d1 = flatten(form["seq1"].value)
    d2 = flatten(form["seq2"].value)

    check1 = clean(d1)
    check2 = clean(d2)

    print "Content-type: text/html"

    if (check1 == False):
        print " CANNOT PROCESS ORIGINAL SEQUENCE "
        print "<p />"
        print " Error 1: Sequence Length not divisible by 3"
        print "<p />"
        return

    elif (check2 == False):
        print " CANNOT PROCESS MUTATED SEQUENCE "
        print "<p />"
        print " Error 1: Sequence Length not divisible by 3"
        print "<p />"
        return

    else:

        data1 = convertDNAtoRNA(d1)
        data2 = convertDNAtoRNA(d2)

        Ks, Ka = KaKsCalcNG(data1, data2)
        Ka2, Ks2, va, vs = KaKsCalcLWL(data1, data2)  # include Variance!!!!
        Ka3, Ks3 = KaKsCalcMLWL(data1, data2)

        values = {}
        values["KaNG"] = Ka
        values["KsNG"] = Ks
        values["KaKsNG"] = Ka/Ks

        values["KaLWL"] = Ka2
        values["KsLWL"] = Ks2
        values["KaKsLWL"] = Ka2/Ks2

        values["KaMLWL"] = Ka3
        values["KsMLWL"] = Ks3
        values["KaKsMLWL"] = Ka3/Ks3

        pageFile = file("results.html")
        pageHtml = pageFile.read()
        pageFile.close()

        print "content-type: text/html"
        print
        print pageHtml % values


main()
