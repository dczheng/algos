#!/usr/bin/env python3
import re
import sys
    
def gen_protos( f ):

    ms = [ '#include', '#define', '#if', '#un', '#endif' ]
    lines = open( f ).readlines()

    ll = []
    for l in lines:
        if l.split() == []:
            continue
        ll.append( l )
    lines = ll

    n = len(lines)
    i = 0
    ss = []
    while i < n:
        '''
        remove macros and `\\` comments
        '''
        l = lines[i].strip()

        flag = 0
        for m in ms:
            if m in l:
                flag = 1

        i += 1
        if flag:
            if l[-1] != '\\':
                continue

            for x in range(i, n):
                t = lines[x].strip()
                if t[-1] == '\\':
                    i += 1
                else:
                    break
            i += 1
            continue

        t = re.search( '//', l )
        if t:
            ss.append( l[:t.span()[0]] )
        else:
            ss.append( l )
    s = ''.join( ss ) 

    def remove_some_char( s, c0, c1, e=None, l=1 ):
        '''remove block'''
        ss = []
        n = len(s)
        i = 0
        while i < n:
            c = s[i:i+l]
            if c == c0 or c == c1:
                i += l 
            else:
                c = s[i]
                i += 1

            if c != c1:
                ss.append( c )
                continue

            cc = ss.pop()
            while cc != c0 and len(ss) != 0:
                cc = ss.pop()
            if e:
               ss.append( e )
        return ''.join(ss)

    s = remove_some_char( s, '{', '}', e=';' )
    s = remove_some_char( s, r'/*', r'*/', l=2  )

    s = s.split( ';' )
    ss = []
    for l in s:
        if '(' in l:
            ss.append( l.strip() + ';' )

    return ss
    

ss = gen_protos( sys.argv[1] )
f = open( sys.argv[1][:-2] + '.h', 'w' )
for l in ss:
    f.write( "%s\n"%l)
f.close()
