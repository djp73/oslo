# Trying to sort out the problems in afsr2

attach afsr.sage

def uktest(q,p,N,t):
    gamma = 1
    for i in range(q):
        if i*p%q==1:
            gamma = i
    rho = 1 
    for i in range(p):
        if i*q%p==1:
            rho = i
    
#    print "Test 1"
#    print adic_seq(u,q,p,N)
#    for i in range(N):
#        print adic_seq(useq[i],q,p,N)

#    print "\nTest 2"
    csum=0
    dsum=0
    for k in range(t):
        u=0
        while (u==q or u==0):
            u=randint(-q,q)
        useq=[(u*gamma**j)%q-q for j in range(N)]
        alpha=adic_seq(u,q,p,N)
#        print alpha
#        aseq = [rho*useq[j]%p for j in range(N)]
        aseq = adic_seq2(u,q,p,N)
        r=ceil(log(q,p))
        c=1
        d=1
        rs=[]
        for i in range(r,N):
            if alpha[2][i]<>aseq[2][i]:
                c=0
        for i in range(r):
            if alpha[2][i]<>aseq[2][i]:
                d=0
            else:
                rs+=[i]
        if c==0:
            print "wrong after r"
            print alpha
            print aseq
        if d==0:
            print "u/q=",u,"/",q," digits correct before r",rs
#           print alpha
#           print aseq
        csum+=c
        dsum+=d
    print "q=",q," p=",p,csum==t,dsum==t
    return

# Results
    """

    sage: for i in range(500):   
    ....:     q=randint(2,2^10);p=randint(2,2^10);
    ....:     if gcd(q,p)==1:
    ....:         uktest(q,p,100,20)
    ....:         
    q= 450  p= 919 True True
    q= 256  p= 883 True True
    q= 217  p= 660 True True
    q= 750  p= 533 True True
    q= 87  p= 713 True True
    q= 247  p= 678 True True
    q= 556  p= 663 True True
    q= 163  p= 519 True True
    q= 447  p= 568 True True
    q= 1014  p= 319 True True
    q= 1005  p= 79 True True
    q= 68  p= 953 True True
    q= 446  p= 699 True True
    q= 458  p= 1013 True True
    q= 509  p= 218 True True
    q= 585  p= 193 True True
    q= 171  p= 251 True True
    q= 626  p= 449 True True
    q= 1018  p= 723 True True
    q= 781  p= 525 True True
    q= 717  p= 89 True True
    q= 991  p= 24 True True
    q= 795  p= 398 True True
    q= 442  p= 159 True True
    q= 862  p= 593 True True
    q= 742  p= 629 True True
    q= 516  p= 115 True True
    q= 13  p= 6 True True
    q= 949  p= 519 True True
    q= 380  p= 499 True True
    q= 535  p= 943 True True
    q= 713  p= 822 True True
    q= 712  p= 701 True True
    q= 440  p= 607 True True
    q= 349  p= 952 True True
    q= 479  p= 913 True True
    q= 516  p= 947 True True
    q= 813  p= 998 True True
    q= 863  p= 596 True True
    q= 694  p= 215 True True
    q= 221  p= 627 True True
    q= 300  p= 341 True True
    q= 3  p= 931 True True
    q= 307  p= 324 True True
    q= 917  p= 698 True True
    q= 89  p= 556 True True
    q= 131  p= 884 True True
    q= 585  p= 797 True True
    q= 841  p= 185 True True
    q= 919  p= 889 True True
    q= 220  p= 571 True True
    q= 799  p= 474 True True
    q= 428  p= 367 True True
    q= 196  p= 243 True True
    q= 599  p= 987 True True
    q= 193  p= 186 True True
    q= 343  p= 320 True True
    q= 148  p= 125 True True
    q= 766  p= 373 True True
    q= 527  p= 991 True True
    q= 17  p= 359 True True
    q= 799  p= 193 True True
    q= 781  p= 998 True True
    q= 443  p= 848 True True
    q= 161  p= 556 True True
    q= 961  p= 743 True True
    q= 670  p= 59 True True
    q= 303  p= 461 True True
    q= 811  p= 542 True True
    q= 293  p= 123 True True
    q= 347  p= 292 True True
    q= 850  p= 403 True True
    q= 73  p= 254 True True
    q= 613  p= 517 True True
    q= 759  p= 629 True True
    q= 99  p= 733 True True
    q= 701  p= 508 True True
    q= 961  p= 769 True True
    q= 435  p= 503 True True
    q= 451  p= 983 True True
    q= 909  p= 896 True True
    q= 626  p= 283 True True
    q= 981  p= 503 True True
    q= 171  p= 854 True True
    q= 563  p= 785 True True
    q= 754  p= 727 True True
    q= 401  p= 831 True True
    q= 775  p= 482 True True
    q= 983  p= 754 True True
    q= 403  p= 711 True True
    q= 224  p= 373 True True
    q= 15  p= 481 True True
    q= 331  p= 963 True True
    q= 849  p= 574 True True
    q= 845  p= 379 True True
    q= 979  p= 1006 True True
    q= 742  p= 307 True True
    q= 116  p= 751 True True
    q= 783  p= 280 True True
    q= 613  p= 99 True True
    q= 799  p= 318 True True
    q= 774  p= 281 True True
    q= 40  p= 771 True True
    q= 657  p= 782 True True
    q= 25  p= 1011 True True
    q= 82  p= 481 True True
    q= 202  p= 71 True True
    q= 855  p= 761 True True
    q= 445  p= 347 True True
    q= 243  p= 232 True True
    q= 323  p= 314 True True
    q= 517  p= 502 True True
    q= 80  p= 841 True True
    q= 457  p= 639 True True
    q= 139  p= 86 True True
    q= 995  p= 829 True True
    q= 815  p= 269 True True
    q= 123  p= 472 True True
    q= 737  p= 750 True True
    q= 319  p= 225 True True
    q= 279  p= 485 True True
    q= 838  p= 611 True True
    q= 283  p= 665 True True
    q= 41  p= 202 True True
    q= 179  p= 332 True True
    q= 35  p= 942 True True
    q= 572  p= 339 True True
    q= 512  p= 917 True True
    q= 158  p= 167 True True
    q= 772  p= 825 True True
    q= 111  p= 412 True True
    q= 223  p= 997 True True
    q= 1004  p= 741 True True
    q= 47  p= 317 True True
    q= 41  p= 517 True True
    q= 79  p= 730 True True
    q= 314  p= 495 True True
    q= 690  p= 953 True True
    q= 367  p= 377 True True
    q= 832  p= 471 True True
    q= 681  p= 800 True True
    q= 814  p= 69 True True
    q= 493  p= 664 True True
    q= 967  p= 262 True True
    q= 997  p= 355 True True
    q= 615  p= 518 True True
    q= 409  p= 173 True True
    q= 891  p= 163 True True
    q= 694  p= 379 True True
    q= 879  p= 614 True True
    q= 886  p= 705 True True
    q= 613  p= 989 True True
    q= 479  p= 17 True True
    q= 961  p= 92 True True
    q= 925  p= 1004 True True
    q= 179  p= 23 True True
    q= 283  p= 705 True True
    q= 203  p= 170 True True
    q= 421  p= 685 True True
    q= 757  p= 309 True True
    q= 137  p= 362 True True
    q= 487  p= 39 True True
    q= 775  p= 427 True True
    q= 314  p= 397 True True
    q= 909  p= 773 True True
    q= 589  p= 327 True True
    q= 326  p= 381 True True
    q= 737  p= 354 True True
    q= 227  p= 605 True True
    q= 589  p= 366 True True
    q= 671  p= 436 True True
    q= 791  p= 243 True True
    q= 1013  p= 512 True True
    q= 181  p= 150 True True
    q= 997  p= 213 True True
    q= 943  p= 840 True True
    q= 408  p= 293 True True
    q= 113  p= 203 True True
    q= 878  p= 675 True True
    q= 869  p= 109 True True
    q= 684  p= 617 True True
    q= 809  p= 733 True True
    q= 539  p= 901 True True
    q= 946  p= 291 True True
    q= 691  p= 923 True True
    q= 797  p= 392 True True
    q= 793  p= 84 True True
    q= 116  p= 269 True True
    q= 908  p= 83 True True
    q= 23  p= 881 True True
    q= 8  p= 207 True True
    q= 559  p= 161 True True
    q= 507  p= 200 True True
    q= 261  p= 22 True True
    q= 618  p= 79 True True
    q= 273  p= 341 True True
    q= 398  p= 359 True True
    q= 857  p= 1006 True True
    q= 842  p= 489 True True
    q= 382  p= 161 True True
    q= 161  p= 755 True True
    q= 975  p= 568 True True
    q= 788  p= 581 True True
    q= 249  p= 575 True True
    q= 392  p= 907 True True
    q= 262  p= 133 True True
    q= 883  p= 1003 True True
    q= 690  p= 941 True True
    q= 347  p= 411 True True
    q= 126  p= 221 True True
    q= 305  p= 406 True True
    q= 837  p= 790 True True
    q= 887  p= 775 True True
    q= 649  p= 24 True True
    q= 971  p= 564 True True
    q= 599  p= 764 True True
    q= 284  p= 535 True True
    q= 499  p= 112 True True
    q= 847  p= 635 True True
    q= 205  p= 499 True True
    q= 272  p= 11 True True
    q= 533  p= 649 True True
    q= 521  p= 535 True True
    q= 50  p= 537 True True
    q= 939  p= 41 True True
    q= 769  p= 624 True True
    q= 195  p= 857 True True
    q= 860  p= 221 True True
    q= 739  p= 675 True True
    q= 677  p= 586 True True
    q= 907  p= 727 True True
    q= 871  p= 599 True True
    q= 728  p= 965 True True
    q= 803  p= 362 True True
    q= 589  p= 839 True True
    q= 351  p= 421 True True
    q= 810  p= 229 True True
    q= 427  p= 299 True True
    q= 797  p= 774 True True
    q= 467  p= 193 True True
    q= 72  p= 487 True True
    q= 394  p= 327 True True
    q= 914  p= 141 True True
    q= 37  p= 298 True True
    q= 808  p= 51 True True
    q= 863  p= 343 True True
    q= 637  p= 493 True True
    q= 315  p= 4 True True
    q= 62  p= 719 True True
    q= 499  p= 401 True True
    q= 671  p= 272 True True
    q= 653  p= 741 True True
    q= 31  p= 792 True True
    q= 164  p= 161 True True
    q= 748  p= 1011 True True
    q= 207  p= 731 True True
    q= 309  p= 212 True True
    q= 622  p= 449 True True
    q= 837  p= 716 True True
    q= 409  p= 1001 True True
    q= 853  p= 327 True True
    q= 513  p= 397 True True
    q= 47  p= 919 True True
    q= 39  p= 1000 True True
    q= 567  p= 370 True True
    q= 229  p= 448 True True
    q= 395  p= 614 True True
    q= 960  p= 661 True True
    q= 994  p= 121 True True
    q= 73  p= 397 True True
    q= 95  p= 293 True True
    q= 484  p= 181 True True
    q= 562  p= 319 True True
    q= 267  p= 388 True True
    q= 871  p= 74 True True
    q= 54  p= 299 True True
    q= 406  p= 701 True True
    q= 539  p= 905 True True
    q= 796  p= 925 True True
    q= 23  p= 197 True True
    q= 219  p= 211 True True
    q= 271  p= 887 True True
    q= 232  p= 807 True True
    q= 269  p= 196 True True
    q= 437  p= 801 True True
    q= 445  p= 891 True True
    q= 208  p= 153 True True
    q= 267  p= 896 True True
    q= 753  p= 979 True True
    q= 441  p= 752 True True
    q= 544  p= 65 True True
    q= 993  p= 4 True True
    q= 193  p= 803 True True
    q= 862  p= 265 True True
    q= 129  p= 676 True True
    q= 487  p= 357 True True
