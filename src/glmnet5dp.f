c     mortran 2.0     (version of 7/04/75 mod 7/4/87 (ajc))              
      subroutine get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)          772
      implicit double precision(a-h,o-z)                                    773
      data sml0,eps0,big0,mnlam0,rsqmax0,pmin0,exmx0  /1.0d-5,1.0d-6,9.9    775 
     *d35,5,0.999,1.0d-9,250.0/
      sml=sml0                                                              775
      eps=eps0                                                              775
      big=big0                                                              775
      mnlam=mnlam0                                                          775
      rsqmax=rsqmax0                                                        776
      pmin=pmin0                                                            776
      exmx=exmx0                                                            777
      return                                                                778
      entry chg_fract_dev(arg)                                              778
      sml0=arg                                                              778
      return                                                                779
      entry chg_dev_max(arg)                                                779
      rsqmax0=arg                                                           779
      return                                                                780
      entry chg_min_flmin(arg)                                              780
      eps0=arg                                                              780
      return                                                                781
      entry chg_big(arg)                                                    781
      big0=arg                                                              781
      return                                                                782
      entry chg_min_lambdas(irg)                                            782
      mnlam0=irg                                                            782
      return                                                                783
      entry chg_min_null_prob(arg)                                          783
      pmin0=arg                                                             783
      return                                                                784
      entry chg_max_exp(arg)                                                784
      exmx0=arg                                                             784
      return                                                                785
      end                                                                   786
      subroutine elnet  (ka,parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,u    789 
     *lam,thr,isd,intr,maxit,  beta0,isg,plam,lmu,a0,ca,ia,nin,rsq,alm,n
     *lp,jerr)
      implicit double precision(a-h,o-z)                                    790
      double precision x(no,ni),y(no),w(no),vp(ni),ca(nx,nlam),cl(2,ni)     791
      double precision beta0(ni)                                            792
      double precision ulam(nlam),a0(nlam),rsq(nlam),alm(nlam)              793
      integer jd(*),ia(nx),nin(nlam)                                        794
      double precision, dimension (:), allocatable :: vq;                       
      if(maxval(vp) .gt. 0.0)goto 10021                                     797
      jerr=10000                                                            797
      return                                                                797
10021 continue                                                              798
      allocate(vq(1:ni),stat=jerr)                                          798
      if(jerr.ne.0) return                                                  799
      vq=max(0d0,vp)                                                        799
      vq=vq*ni/sum(vq)                                                      800
      if(ka .ne. 1)goto 10041                                               801
      call elnetu  (parm,no,ni,x,y,w,jd,vq,cl,ne,nx,nlam,flmin,ulam,thr,    804 
     *isd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 10051                                                            805
10041 continue                                                              806
      call elnetn (parm,no,ni,x,y,w,jd,vq,cl,ne,nx,nlam,flmin,ulam,thr,i    809 
     *sd,intr,maxit,  beta0,isg,plam,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
10051 continue                                                              810
10031 continue                                                              810
      deallocate(vq)                                                        811
      return                                                                812
      end                                                                   813
      subroutine elnetu  (parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,ula    816 
     *m,thr,isd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                    817
      double precision x(no,ni),y(no),w(no),vp(ni),ulam(nlam),cl(2,ni)      818
      double precision ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)             819
      integer jd(*),ia(nx),nin(nlam)                                        820
      double precision, dimension (:), allocatable :: xm,xs,g,xv,vlam           
      integer, dimension (:), allocatable :: ju                                 
      allocate(g(1:ni),stat=jerr)                                           825
      if(jerr.ne.0) return                                                  826
      allocate(xm(1:ni),stat=jerr)                                          827
      if(jerr.ne.0) return                                                  828
      allocate(xs(1:ni),stat=jerr)                                          829
      if(jerr.ne.0) return                                                  830
      allocate(ju(1:ni),stat=jerr)                                          831
      if(jerr.ne.0) return                                                  832
      allocate(xv(1:ni),stat=jerr)                                          833
      if(jerr.ne.0) return                                                  834
      allocate(vlam(1:nlam),stat=jerr)                                      835
      if(jerr.ne.0) return                                                  836
      call chkvars(no,ni,x,ju)                                              837
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  838
      if(maxval(ju) .gt. 0)goto 10071                                       838
      jerr=7777                                                             838
      return                                                                838
10071 continue                                                              839
      call standard(no,ni,x,y,w,isd,intr,ju,g,xm,xs,ym,ys,xv,jerr)          840
      if(jerr.ne.0) return                                                  841
      cl=cl/ys                                                              841
      if(isd .le. 0)goto 10091                                              841
10100 do 10101 j=1,ni                                                       841
      cl(:,j)=cl(:,j)*xs(j)                                                 841
10101 continue                                                              841
10102 continue                                                              841
10091 continue                                                              842
      if(flmin.ge.1.0) vlam=ulam/ys                                         843
      call elnet1(parm,ni,ju,vp,cl,g,no,ne,nx,x,nlam,flmin,vlam,thr,maxi    845 
     *t,xv,  lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                  846
10110 do 10111 k=1,lmu                                                      846
      alm(k)=ys*alm(k)                                                      846
      nk=nin(k)                                                             847
10120 do 10121 l=1,nk                                                       847
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                          847
10121 continue                                                              847
10122 continue                                                              847
      a0(k)=0.0                                                             848
      if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))           849
10111 continue                                                              850
10112 continue                                                              850
      deallocate(xm,xs,g,ju,xv,vlam)                                        851
      return                                                                852
      end                                                                   853
      subroutine standard (no,ni,x,y,w,isd,intr,ju,g,xm,xs,ym,ys,xv,jerr    854 
     *)
      implicit double precision(a-h,o-z)                                    855
      double precision x(no,ni),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)      856
      integer ju(ni)                                                        857
      double precision, dimension (:), allocatable :: v                         
      allocate(v(1:no),stat=jerr)                                           860
      if(jerr.ne.0) return                                                  861
      w=w/sum(w)                                                            861
      v=sqrt(w)                                                             862
      if(intr .ne. 0)goto 10141                                             862
      ym=0.0                                                                862
      y=v*y                                                                 863
      ys=sqrt(dot_product(y,y)-dot_product(v,y)**2)                         863
      y=y/ys                                                                864
10150 do 10151 j=1,ni                                                       864
      if(ju(j).eq.0)goto 10151                                              864
      xm(j)=0.0                                                             864
      x(:,j)=v*x(:,j)                                                       865
      xv(j)=dot_product(x(:,j),x(:,j))                                      866
      if(isd .eq. 0)goto 10171                                              866
      xbq=dot_product(v,x(:,j))**2                                          866
      vc=xv(j)-xbq                                                          867
      xs(j)=sqrt(vc)                                                        867
      x(:,j)=x(:,j)/xs(j)                                                   867
      xv(j)=1.0+xbq/vc                                                      868
      goto 10181                                                            869
10171 continue                                                              869
      xs(j)=1.0                                                             869
10181 continue                                                              870
10161 continue                                                              870
10151 continue                                                              871
10152 continue                                                              871
      goto 10191                                                            872
10141 continue                                                              873
10200 do 10201 j=1,ni                                                       873
      if(ju(j).eq.0)goto 10201                                              874
      xm(j)=dot_product(w,x(:,j))                                           874
      x(:,j)=v*(x(:,j)-xm(j))                                               875
      xv(j)=dot_product(x(:,j),x(:,j))                                      875
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                        876
10201 continue                                                              877
10202 continue                                                              877
      if(isd .ne. 0)goto 10221                                              877
      xs=1.0                                                                877
      goto 10231                                                            878
10221 continue                                                              879
10240 do 10241 j=1,ni                                                       879
      if(ju(j).eq.0)goto 10241                                              879
      x(:,j)=x(:,j)/xs(j)                                                   879
10241 continue                                                              880
10242 continue                                                              880
      xv=1.0                                                                881
10231 continue                                                              882
10211 continue                                                              882
      ym=dot_product(w,y)                                                   882
      y=v*(y-ym)                                                            882
      ys=sqrt(dot_product(y,y))                                             882
      y=y/ys                                                                883
10191 continue                                                              884
10131 continue                                                              884
      g=0.0                                                                 884
10250 do 10251 j=1,ni                                                       884
      if(ju(j).ne.0) g(j)=dot_product(y,x(:,j))                             884
10251 continue                                                              885
10252 continue                                                              885
      deallocate(v)                                                         886
      return                                                                887
      end                                                                   888
      subroutine elnet1 (beta,ni,ju,vp,cl,g,no,ne,nx,x,nlam,flmin,ulam,t    890 
     *hr,maxit,xv,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      implicit double precision(a-h,o-z)                                    891
      double precision vp(ni),g(ni),x(no,ni),ulam(nlam),ao(nx,nlam)         892
      double precision rsqo(nlam),almo(nlam),xv(ni)                         893
      double precision cl(2,ni)                                             894
      integer ju(ni),ia(nx),kin(nlam)                                       895
      double precision, dimension (:), allocatable :: a,da                      
      integer, dimension (:), allocatable :: mm                                 
      double precision, dimension (:,:), allocatable :: c                       
      allocate(c(1:ni,1:nx),stat=jerr)                                          
      if(jerr.ne.0) return;                                                     
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)                903
      allocate(a(1:ni),stat=jerr)                                           904
      if(jerr.ne.0) return                                                  905
      allocate(mm(1:ni),stat=jerr)                                          906
      if(jerr.ne.0) return                                                  907
      allocate(da(1:ni),stat=jerr)                                          908
      if(jerr.ne.0) return                                                  909
      bta=beta                                                              909
      omb=1.0-bta                                                           910
      if(flmin .ge. 1.0)goto 10271                                          910
      eqs=max(eps,flmin)                                                    910
      alf=eqs**(1.0/(nlam-1))                                               910
10271 continue                                                              911
      rsq=0.0                                                               911
      a=0.0                                                                 911
      mm=0                                                                  911
      nlp=0                                                                 911
      nin=nlp                                                               911
      iz=0                                                                  911
      mnl=min(mnlam,nlam)                                                   913
      alm=0.0                                                               915
10280 do 10281 m=1,nlam                                                     916
      if(flmin .lt. 1.0)goto 10301                                          916
      alm=ulam(m)                                                           916
      goto 10291                                                            917
10301 if(m .le. 2)goto 10311                                                917
      alm=alm*alf                                                           917
      goto 10291                                                            918
10311 if(m .ne. 1)goto 10321                                                918
      alm=big                                                               918
      goto 10331                                                            919
10321 continue                                                              919
      alm=0.0                                                               920
10340 do 10341 j=1,ni                                                       920
      if(ju(j).eq.0)goto 10341                                              920
      if(vp(j).le.0.0)goto 10341                                            921
      alm=max(alm,abs(g(j))/vp(j))                                          922
10341 continue                                                              923
10342 continue                                                              923
      alm=alf*alm/max(bta,1.0d-3)                                           924
10331 continue                                                              925
10291 continue                                                              925
      dem=alm*omb                                                           925
      ab=alm*bta                                                            925
      rsq0=rsq                                                              925
      jz=1                                                                  926
10350 continue                                                              926
10351 continue                                                              926
      if(iz*jz.ne.0) go to 10360                                            926
      nlp=nlp+1                                                             926
      dlx=0.0                                                               927
10370 do 10371 k=1,ni                                                       927
      if(ju(k).eq.0)goto 10371                                              928
      ak=a(k)                                                               928
      u=g(k)+ak*xv(k)                                                       928
      v=abs(u)-vp(k)*ab                                                     928
      a(k)=0.0                                                              930
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d    931 
     *em)))
      if(a(k).eq.ak)goto 10371                                              932
      if(mm(k) .ne. 0)goto 10391                                            932
      nin=nin+1                                                             932
      if(nin.gt.nx)goto 10372                                               933
10400 do 10401 j=1,ni                                                       933
      if(ju(j).eq.0)goto 10401                                              934
      if(mm(j) .eq. 0)goto 10421                                            934
      c(j,nin)=c(k,mm(j))                                                   934
      goto 10401                                                            934
10421 continue                                                              935
      if(j .ne. k)goto 10441                                                935
      c(j,nin)=xv(j)                                                        935
      goto 10401                                                            935
10441 continue                                                              936
      c(j,nin)=dot_product(x(:,j),x(:,k))                                   937
10401 continue                                                              938
10402 continue                                                              938
      mm(k)=nin                                                             938
      ia(nin)=k                                                             939
10391 continue                                                              940
      del=a(k)-ak                                                           940
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                      941
      dlx=max(xv(k)*del**2,dlx)                                             942
10450 do 10451 j=1,ni                                                       942
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del                               942
10451 continue                                                              943
10452 continue                                                              943
10371 continue                                                              944
10372 continue                                                              944
      if(dlx.lt.thr)goto 10352                                              944
      if(nin.gt.nx)goto 10352                                               945
      if(nlp .le. maxit)goto 10471                                          945
      jerr=-m                                                               945
      return                                                                945
10471 continue                                                              946
10360 continue                                                              946
      iz=1                                                                  946
      da(1:nin)=a(ia(1:nin))                                                947
10480 continue                                                              947
10481 continue                                                              947
      nlp=nlp+1                                                             947
      dlx=0.0                                                               948
10490 do 10491 l=1,nin                                                      948
      k=ia(l)                                                               948
      ak=a(k)                                                               948
      u=g(k)+ak*xv(k)                                                       948
      v=abs(u)-vp(k)*ab                                                     949
      a(k)=0.0                                                              951
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d    952 
     *em)))
      if(a(k).eq.ak)goto 10491                                              953
      del=a(k)-ak                                                           953
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                      954
      dlx=max(xv(k)*del**2,dlx)                                             955
10500 do 10501 j=1,nin                                                      955
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                                  955
10501 continue                                                              956
10502 continue                                                              956
10491 continue                                                              957
10492 continue                                                              957
      if(dlx.lt.thr)goto 10482                                              957
      if(nlp .le. maxit)goto 10521                                          957
      jerr=-m                                                               957
      return                                                                957
10521 continue                                                              958
      goto 10481                                                            959
10482 continue                                                              959
      da(1:nin)=a(ia(1:nin))-da(1:nin)                                      960
10530 do 10531 j=1,ni                                                       960
      if(mm(j).ne.0)goto 10531                                              961
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))            962
10531 continue                                                              963
10532 continue                                                              963
      jz=0                                                                  964
      goto 10351                                                            965
10352 continue                                                              965
      if(nin .le. nx)goto 10551                                             965
      jerr=-10000-m                                                         965
      goto 10282                                                            965
10551 continue                                                              966
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                 966
      kin(m)=nin                                                            967
      rsqo(m)=rsq                                                           967
      almo(m)=alm                                                           967
      lmu=m                                                                 968
      if(m.lt.mnl)goto 10281                                                968
      if(flmin.ge.1.0)goto 10281                                            969
      me=0                                                                  969
10560 do 10561 j=1,nin                                                      969
      if(ao(j,m).ne.0.0) me=me+1                                            969
10561 continue                                                              969
10562 continue                                                              969
      if(me.gt.ne)goto 10282                                                970
      if(rsq-rsq0.lt.sml*rsq)goto 10282                                     970
      if(rsq.gt.rsqmax)goto 10282                                           971
10281 continue                                                              972
10282 continue                                                              972
      deallocate(a,mm,c,da)                                                 973
      return                                                                974
      end                                                                   975
      subroutine elnetn (parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam    977 
     *,thr,isd,  intr,maxit,beta0,isg,plam,lmu,a0,ca,ia,nin,rsq,alm,nlp,
     *jerr)
      implicit double precision(a-h,o-z)                                    978
      double precision vp(ni),x(no,ni),y(no),w(no),ulam(nlam),cl(2,ni)      979
      double precision beta0(ni)                                            980
      double precision ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)             981
      integer jd(*),ia(nx),nin(nlam)                                        982
      double precision, dimension (:), allocatable :: xm,xs,xv,vlam             
      integer, dimension (:), allocatable :: ju                                 
      allocate(xm(1:ni),stat=jerr)                                          987
      if(jerr.ne.0) return                                                  988
      allocate(xs(1:ni),stat=jerr)                                          989
      if(jerr.ne.0) return                                                  990
      allocate(ju(1:ni),stat=jerr)                                          991
      if(jerr.ne.0) return                                                  992
      allocate(xv(1:ni),stat=jerr)                                          993
      if(jerr.ne.0) return                                                  994
      allocate(vlam(1:nlam),stat=jerr)                                      995
      if(jerr.ne.0) return                                                  996
      call chkvars(no,ni,x,ju)                                              997
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  998
      if(maxval(ju) .gt. 0)goto 10581                                       998
      jerr=7777                                                             998
      return                                                                998
10581 continue                                                              999
      call standard1(no,ni,x,y,w,isd,intr,ju,xm,xs,ym,ys,xv,jerr)          1000
      if(jerr.ne.0) return                                                 1001
      beta0=beta0/ys                                                       1002
      cl=cl/ys                                                             1002
      if(isd .le. 0)goto 10601                                             1002
10610 do 10611 j=1,ni                                                      1002
      cl(:,j)=cl(:,j)*xs(j)                                                1003
      beta0(j)=beta0(j)*xs(j)                                              1003
10611 continue                                                             1003
10612 continue                                                             1003
10601 continue                                                             1004
      if(flmin .lt. 1.0)goto 10631                                         1004
      vlam=ulam/ys                                                         1004
      plam=plam/ys                                                         1004
10631 continue                                                             1005
      call elnet2(parm,ni,ju,vp,cl,y,no,ne,nx,x,nlam,flmin,vlam,thr,maxi   1007 
     *t,  beta0,isg,plam,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 1008
10640 do 10641 k=1,lmu                                                     1008
      alm(k)=ys*alm(k)                                                     1008
      nk=nin(k)                                                            1009
10650 do 10651 l=1,nk                                                      1009
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                         1009
10651 continue                                                             1009
10652 continue                                                             1009
      a0(k)=0.0                                                            1010
      if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))          1011
10641 continue                                                             1012
10642 continue                                                             1012
      deallocate(xm,xs,ju,xv,vlam)                                         1013
      return                                                               1014
      end                                                                  1015
      subroutine standard1 (no,ni,x,y,w,isd,intr,ju,xm,xs,ym,ys,xv,jerr)   1016
      implicit double precision(a-h,o-z)                                   1017
      double precision x(no,ni),y(no),w(no),xm(ni),xs(ni),xv(ni)           1017
      integer ju(ni)                                                       1018
      double precision, dimension (:), allocatable :: v                         
      allocate(v(1:no),stat=jerr)                                          1021
      if(jerr.ne.0) return                                                 1022
      w=w/sum(w)                                                           1022
      v=sqrt(w)                                                            1023
      if(intr .ne. 0)goto 10671                                            1023
      ym=0.0                                                               1023
      y=v*y                                                                1024
      ys=sqrt(dot_product(y,y)-dot_product(v,y)**2)                        1024
      y=y/ys                                                               1025
10680 do 10681 j=1,ni                                                      1025
      if(ju(j).eq.0)goto 10681                                             1025
      xm(j)=0.0                                                            1025
      x(:,j)=v*x(:,j)                                                      1026
      xv(j)=dot_product(x(:,j),x(:,j))                                     1027
      if(isd .eq. 0)goto 10701                                             1027
      xbq=dot_product(v,x(:,j))**2                                         1027
      vc=xv(j)-xbq                                                         1028
      xs(j)=sqrt(vc)                                                       1028
      x(:,j)=x(:,j)/xs(j)                                                  1028
      xv(j)=1.0+xbq/vc                                                     1029
      goto 10711                                                           1030
10701 continue                                                             1030
      xs(j)=1.0                                                            1030
10711 continue                                                             1031
10691 continue                                                             1031
10681 continue                                                             1032
10682 continue                                                             1032
      go to 10720                                                          1033
10671 continue                                                             1034
10730 do 10731 j=1,ni                                                      1034
      if(ju(j).eq.0)goto 10731                                             1035
      xm(j)=dot_product(w,x(:,j))                                          1035
      x(:,j)=v*(x(:,j)-xm(j))                                              1036
      xv(j)=dot_product(x(:,j),x(:,j))                                     1036
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       1037
10731 continue                                                             1038
10732 continue                                                             1038
      if(isd .ne. 0)goto 10751                                             1038
      xs=1.0                                                               1038
      goto 10761                                                           1039
10751 continue                                                             1039
10770 do 10771 j=1,ni                                                      1039
      if(ju(j).eq.0)goto 10771                                             1039
      x(:,j)=x(:,j)/xs(j)                                                  1039
10771 continue                                                             1040
10772 continue                                                             1040
      xv=1.0                                                               1041
10761 continue                                                             1042
10741 continue                                                             1042
      ym=dot_product(w,y)                                                  1042
      y=v*(y-ym)                                                           1042
      ys=sqrt(dot_product(y,y))                                            1042
      y=y/ys                                                               1043
10720 continue                                                             1043
      deallocate(v)                                                        1044
      return                                                               1045
      end                                                                  1046
      subroutine elnet2(beta,ni,ju,vp,cl,y,no,ne,nx,x,nlam,flmin,ulam,th   1048 
     *r,maxit,  beta0,isg,plam,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      implicit double precision(a-h,o-z)                                   1049
      double precision vp(ni),y(no),x(no,ni),ulam(nlam),ao(nx,nlam)        1050
      double precision rsqo(nlam),almo(nlam),xv(ni)                        1051
      double precision cl(2,ni),beta0(ni)                                  1052
      integer ju(ni),ia(nx),kin(nlam)                                      1053
      double precision, dimension (:), allocatable :: a,g                       
      integer, dimension (:), allocatable :: mm,ix                              
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)               1058
      allocate(a(1:ni),stat=jerr)                                          1059
      if(jerr.ne.0) return                                                 1060
      allocate(mm(1:ni),stat=jerr)                                         1061
      if(jerr.ne.0) return                                                 1062
      allocate(g(1:ni),stat=jerr)                                          1063
      if(jerr.ne.0) return                                                 1064
      allocate(ix(1:ni),stat=jerr)                                         1065
      if(jerr.ne.0) return                                                 1066
      bta=beta                                                             1066
      omb=1.0-bta                                                          1066
      ix=0                                                                 1068
      alf=1.0                                                              1070
      if(flmin .ge. 1.0)goto 10791                                         1070
      eqs=max(eps,flmin)                                                   1070
      alf=eqs**(1.0/(nlam-1))                                              1070
10791 continue                                                             1071
      rsq=0.0                                                              1071
      a=0.0                                                                1071
      mm=0                                                                 1071
      nlp=0                                                                1071
      nin=nlp                                                              1071
      iz=0                                                                 1071
      mnl=min(mnlam,nlam)                                                  1071
      alm=0.0                                                              1072
      if(flmin .lt. 1.0)goto 10811                                         1072
10820 do 10821 j=1,ni                                                      1072
      if(ju(j).eq.0)goto 10821                                             1072
      y=y-beta0(j)*x(:,j)                                                  1072
10821 continue                                                             1072
10822 continue                                                             1072
10811 continue                                                             1073
10830 do 10831 j=1,ni                                                      1073
      if(ju(j).eq.0)goto 10831                                             1073
      g(j)=abs(dot_product(y,x(:,j)))                                      1073
10831 continue                                                             1074
10832 continue                                                             1074
10840 do 10841 m=1,nlam                                                    1074
      alm0=alm                                                             1075
      if(flmin .lt. 1.0)goto 10861                                         1075
      alm=ulam(m)                                                          1075
      goto 10851                                                           1076
10861 if(m .le. 2)goto 10871                                               1076
      alm=alm*alf                                                          1076
      goto 10851                                                           1077
10871 if(m .ne. 1)goto 10881                                               1077
      alm=big                                                              1077
      goto 10891                                                           1078
10881 continue                                                             1078
      alm0=0.0                                                             1079
10900 do 10901 j=1,ni                                                      1079
      if(ju(j).eq.0)goto 10901                                             1079
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           1079
10901 continue                                                             1080
10902 continue                                                             1080
      alm0=alm0/max(bta,1.0d-3)                                            1080
      alm=alf*alm0                                                         1081
10891 continue                                                             1082
10851 continue                                                             1082
      dem=alm*omb                                                          1082
      ab=alm*bta                                                           1082
      rsq0=rsq                                                             1082
      jz=1                                                                 1083
      if(flmin .lt. 1.0 .or. m .ne. 1)goto 10921                           1084
10930 do 10931 k=1,ni                                                      1084
      if(ju(k).eq.0)goto 10931                                             1084
      a(k)=beta0(k)                                                        1085
      if(abs(beta0(k)) .le. 1e-12)goto 10951                               1085
      ix(k)=1                                                              1085
      nin=nin+1                                                            1085
      mm(k)=nin                                                            1085
      ia(nin)=k                                                            1085
10951 continue                                                             1086
10931 continue                                                             1087
10932 continue                                                             1087
      if(isg .ne. 1)goto 10971                                             1087
      tlam=bta*(2.0*alm-plam)                                              1088
10980 do 10981 k=1,ni                                                      1088
      if(ix(k).eq.1)goto 10981                                             1088
      if(ju(k).eq.0)goto 10981                                             1089
      if(g(k).gt.tlam*vp(k)) ix(k)=1                                       1090
10981 continue                                                             1091
10982 continue                                                             1091
10971 continue                                                             1092
      goto 10991                                                           1093
10921 continue                                                             1094
      tlam=bta*(2.0*alm-alm0)                                              1095
11000 do 11001 k=1,ni                                                      1095
      if(ix(k).eq.1)goto 11001                                             1095
      if(ju(k).eq.0)goto 11001                                             1096
      if(g(k).gt.tlam*vp(k)) ix(k)=1                                       1097
11001 continue                                                             1098
11002 continue                                                             1098
10991 continue                                                             1099
10911 continue                                                             1099
11010 continue                                                             1099
11011 continue                                                             1099
      if(iz*jz.ne.0) go to 10360                                           1100
11020 continue                                                             1100
      nlp=nlp+1                                                            1100
      dlx=0.0                                                              1101
11030 do 11031 k=1,ni                                                      1101
      if(ix(k).eq.0)goto 11031                                             1101
      gk=dot_product(y,x(:,k))                                             1102
      ak=a(k)                                                              1102
      u=gk+ak*xv(k)                                                        1102
      v=abs(u)-vp(k)*ab                                                    1102
      a(k)=0.0                                                             1104
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   1105 
     *em)))
      if(a(k).eq.ak)goto 11031                                             1106
      if(mm(k) .ne. 0)goto 11051                                           1106
      nin=nin+1                                                            1106
      if(nin.gt.nx)goto 11032                                              1107
      mm(k)=nin                                                            1107
      ia(nin)=k                                                            1108
11051 continue                                                             1109
      del=a(k)-ak                                                          1109
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1110
      y=y-del*x(:,k)                                                       1110
      dlx=max(xv(k)*del**2,dlx)                                            1111
11031 continue                                                             1112
11032 continue                                                             1112
      if(nin.gt.nx)goto 11012                                              1113
      if(dlx .ge. thr)goto 11071                                           1113
      ixx=0                                                                1114
11080 do 11081 k=1,ni                                                      1114
      if(ix(k).eq.1)goto 11081                                             1114
      if(ju(k).eq.0)goto 11081                                             1115
      g(k)=abs(dot_product(y,x(:,k)))                                      1116
      if(g(k) .le. ab*vp(k))goto 11101                                     1116
      ix(k)=1                                                              1116
      ixx=1                                                                1116
11101 continue                                                             1117
11081 continue                                                             1118
11082 continue                                                             1118
      if(ixx.eq.1) go to 11020                                             1119
      goto 11012                                                           1120
11071 continue                                                             1121
      if(nlp .le. maxit)goto 11121                                         1121
      jerr=-m                                                              1121
      return                                                               1121
11121 continue                                                             1122
10360 continue                                                             1122
      iz=1                                                                 1123
11130 continue                                                             1123
11131 continue                                                             1123
      nlp=nlp+1                                                            1123
      dlx=0.0                                                              1124
11140 do 11141 l=1,nin                                                     1124
      k=ia(l)                                                              1124
      gk=dot_product(y,x(:,k))                                             1125
      ak=a(k)                                                              1125
      u=gk+ak*xv(k)                                                        1125
      v=abs(u)-vp(k)*ab                                                    1125
      a(k)=0.0                                                             1127
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   1128 
     *em)))
      if(a(k).eq.ak)goto 11141                                             1129
      del=a(k)-ak                                                          1129
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1130
      y=y-del*x(:,k)                                                       1130
      dlx=max(xv(k)*del**2,dlx)                                            1131
11141 continue                                                             1132
11142 continue                                                             1132
      if(dlx.lt.thr)goto 11132                                             1132
      if(nlp .le. maxit)goto 11161                                         1132
      jerr=-m                                                              1132
      return                                                               1132
11161 continue                                                             1133
      goto 11131                                                           1134
11132 continue                                                             1134
      jz=0                                                                 1135
      goto 11011                                                           1136
11012 continue                                                             1136
      if(nin .le. nx)goto 11181                                            1136
      jerr=-10000-m                                                        1136
      goto 10842                                                           1136
11181 continue                                                             1137
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                1137
      kin(m)=nin                                                           1138
      rsqo(m)=rsq                                                          1138
      almo(m)=alm                                                          1138
      lmu=m                                                                1139
      if(m.lt.mnl)goto 10841                                               1139
      if(flmin.ge.1.0)goto 10841                                           1140
      me=0                                                                 1140
11190 do 11191 j=1,nin                                                     1140
      if(ao(j,m).ne.0.0) me=me+1                                           1140
11191 continue                                                             1140
11192 continue                                                             1140
      if(me.gt.ne)goto 10842                                               1141
      if(rsq-rsq0.lt.sml*rsq)goto 10842                                    1141
      if(rsq.gt.rsqmax)goto 10842                                          1142
10841 continue                                                             1143
10842 continue                                                             1143
      deallocate(a,mm,g,ix)                                                1144
      return                                                               1145
      end                                                                  1146
      subroutine chkvars(no,ni,x,ju)                                       1147
      implicit double precision(a-h,o-z)                                   1148
      double precision x(no,ni)                                            1148
      integer ju(ni)                                                       1149
11200 do 11201 j=1,ni                                                      1149
      ju(j)=0                                                              1149
      t=x(1,j)                                                             1150
11210 do 11211 i=2,no                                                      1150
      if(x(i,j).eq.t)goto 11211                                            1150
      ju(j)=1                                                              1150
      goto 11212                                                           1150
11211 continue                                                             1151
11212 continue                                                             1151
11201 continue                                                             1152
11202 continue                                                             1152
      return                                                               1153
      end                                                                  1154
      subroutine uncomp(ni,ca,ia,nin,a)                                    1155
      implicit double precision(a-h,o-z)                                   1156
      double precision ca(*),a(ni)                                         1156
      integer ia(*)                                                        1157
      a=0.0                                                                1157
      if(nin.gt.0) a(ia(1:nin))=ca(1:nin)                                  1158
      return                                                               1159
      end                                                                  1160
      subroutine modval(a0,ca,ia,nin,n,x,f)                                1161
      implicit double precision(a-h,o-z)                                   1162
      double precision ca(nin),x(n,*),f(n)                                 1162
      integer ia(nin)                                                      1163
      f=a0                                                                 1163
      if(nin.le.0) return                                                  1164
11220 do 11221 i=1,n                                                       1164
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                      1164
11221 continue                                                             1165
11222 continue                                                             1165
      return                                                               1166
      end                                                                  1167
      subroutine spelnet  (ka,parm,no,ni,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam   1170 
     *,flmin,ulam,thr,isd,intr,  maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr
     *)
      implicit double precision(a-h,o-z)                                   1171
      double precision x(*),y(no),w(no),vp(ni),ulam(nlam),cl(2,ni)         1172
      double precision ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)            1173
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1174
      double precision, dimension (:), allocatable :: vq;                       
      if(maxval(vp) .gt. 0.0)goto 11241                                    1177
      jerr=10000                                                           1177
      return                                                               1177
11241 continue                                                             1178
      allocate(vq(1:ni),stat=jerr)                                         1178
      if(jerr.ne.0) return                                                 1179
      vq=max(0d0,vp)                                                       1179
      vq=vq*ni/sum(vq)                                                     1180
      if(ka .ne. 1)goto 11261                                              1181
      call spelnetu  (parm,no,ni,x,ix,jx,y,w,jd,vq,cl,ne,nx,nlam,flmin,u   1184 
     *lam,thr,isd,  intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 11271                                                           1185
11261 continue                                                             1186
      call spelnetn (parm,no,ni,x,ix,jx,y,w,jd,vq,cl,ne,nx,nlam,flmin,ul   1189 
     *am,thr,isd,intr,  maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
11271 continue                                                             1190
11251 continue                                                             1190
      deallocate(vq)                                                       1191
      return                                                               1192
      end                                                                  1193
      subroutine spelnetu  (parm,no,ni,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam,f   1196 
     *lmin,ulam,thr,isd,intr,  maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   1197
      double precision x(*),y(no),w(no),vp(ni),ulam(nlam),cl(2,ni)         1198
      double precision ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)            1199
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1200
      double precision, dimension (:), allocatable :: xm,xs,g,xv,vlam           
      integer, dimension (:), allocatable :: ju                                 
      allocate(g(1:ni),stat=jerr)                                          1205
      if(jerr.ne.0) return                                                 1206
      allocate(xm(1:ni),stat=jerr)                                         1207
      if(jerr.ne.0) return                                                 1208
      allocate(xs(1:ni),stat=jerr)                                         1209
      if(jerr.ne.0) return                                                 1210
      allocate(ju(1:ni),stat=jerr)                                         1211
      if(jerr.ne.0) return                                                 1212
      allocate(xv(1:ni),stat=jerr)                                         1213
      if(jerr.ne.0) return                                                 1214
      allocate(vlam(1:nlam),stat=jerr)                                     1215
      if(jerr.ne.0) return                                                 1216
      call spchkvars(no,ni,x,ix,ju)                                        1217
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1218
      if(maxval(ju) .gt. 0)goto 11291                                      1218
      jerr=7777                                                            1218
      return                                                               1218
11291 continue                                                             1219
      call spstandard(no,ni,x,ix,jx,y,w,ju,isd,intr,g,xm,xs,ym,ys,xv,jer   1220 
     *r)
      if(jerr.ne.0) return                                                 1221
      cl=cl/ys                                                             1221
      if(isd .le. 0)goto 11311                                             1221
11320 do 11321 j=1,ni                                                      1221
      cl(:,j)=cl(:,j)*xs(j)                                                1221
11321 continue                                                             1221
11322 continue                                                             1221
11311 continue                                                             1222
      if(flmin.ge.1.0) vlam=ulam/ys                                        1223
      call spelnet1(parm,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,cl,nlam,flmin,vla   1225 
     *m,thr,maxit,  xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 1226
11330 do 11331 k=1,lmu                                                     1226
      alm(k)=ys*alm(k)                                                     1226
      nk=nin(k)                                                            1227
11340 do 11341 l=1,nk                                                      1227
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                         1227
11341 continue                                                             1227
11342 continue                                                             1227
      a0(k)=0.0                                                            1228
      if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))          1229
11331 continue                                                             1230
11332 continue                                                             1230
      deallocate(xm,xs,g,ju,xv,vlam)                                       1231
      return                                                               1232
      end                                                                  1233
      subroutine spstandard (no,ni,x,ix,jx,y,w,ju,isd,intr,g,xm,xs,ym,ys   1234 
     *,xv,jerr)
      implicit double precision(a-h,o-z)                                   1235
      double precision x(*),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)         1236
      integer ix(*),jx(*),ju(ni)                                           1237
      w=w/sum(w)                                                           1238
      if(intr .ne. 0)goto 11361                                            1238
      ym=0.0                                                               1239
      ys=sqrt(dot_product(w,y**2)-dot_product(w,y)**2)                     1239
      y=y/ys                                                               1240
11370 do 11371 j=1,ni                                                      1240
      if(ju(j).eq.0)goto 11371                                             1240
      xm(j)=0.0                                                            1240
      jb=ix(j)                                                             1240
      je=ix(j+1)-1                                                         1241
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)                          1242
      if(isd .eq. 0)goto 11391                                             1242
      xbq=dot_product(w(jx(jb:je)),x(jb:je))**2                            1242
      vc=xv(j)-xbq                                                         1243
      xs(j)=sqrt(vc)                                                       1243
      xv(j)=1.0+xbq/vc                                                     1244
      goto 11401                                                           1245
11391 continue                                                             1245
      xs(j)=1.0                                                            1245
11401 continue                                                             1246
11381 continue                                                             1246
11371 continue                                                             1247
11372 continue                                                             1247
      goto 11411                                                           1248
11361 continue                                                             1249
11420 do 11421 j=1,ni                                                      1249
      if(ju(j).eq.0)goto 11421                                             1250
      jb=ix(j)                                                             1250
      je=ix(j+1)-1                                                         1250
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1251
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 1252
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       1253
11421 continue                                                             1254
11422 continue                                                             1254
      if(isd .ne. 0)goto 11441                                             1254
      xs=1.0                                                               1254
      goto 11451                                                           1254
11441 continue                                                             1254
      xv=1.0                                                               1254
11451 continue                                                             1255
11431 continue                                                             1255
      ym=dot_product(w,y)                                                  1255
      y=y-ym                                                               1255
      ys=sqrt(dot_product(w,y**2))                                         1255
      y=y/ys                                                               1256
11411 continue                                                             1257
11351 continue                                                             1257
      g=0.0                                                                1258
11460 do 11461 j=1,ni                                                      1258
      if(ju(j).eq.0)goto 11461                                             1258
      jb=ix(j)                                                             1258
      je=ix(j+1)-1                                                         1259
      g(j)=dot_product(w(jx(jb:je))*y(jx(jb:je)),x(jb:je))/xs(j)           1260
11461 continue                                                             1261
11462 continue                                                             1261
      return                                                               1262
      end                                                                  1263
      subroutine spelnet1(beta,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,cl,nlam,flm   1265 
     *in,ulam,  thr,maxit,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      implicit double precision(a-h,o-z)                                   1266
      double precision g(ni),vp(ni),x(*),ulam(nlam),w(no)                  1267
      double precision ao(nx,nlam),rsqo(nlam),almo(nlam)                   1268
      double precision xm(ni),xs(ni),xv(ni),cl(2,ni)                       1269
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                          1270
      double precision, dimension (:), allocatable :: a,da                      
      integer, dimension (:), allocatable :: mm                                 
      double precision, dimension (:,:), allocatable :: c                       
      allocate(c(1:ni,1:nx),stat=jerr)                                          
      if(jerr.ne.0) return;                                                     
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)               1278
      allocate(a(1:ni),stat=jerr)                                          1279
      if(jerr.ne.0) return                                                 1280
      allocate(mm(1:ni),stat=jerr)                                         1281
      if(jerr.ne.0) return                                                 1282
      allocate(da(1:ni),stat=jerr)                                         1283
      if(jerr.ne.0) return                                                 1284
      bta=beta                                                             1284
      omb=1.0-bta                                                          1286
      alm=0.0                                                              1286
      alf=1.0                                                              1288
      if(flmin .ge. 1.0)goto 11481                                         1288
      eqs=max(eps,flmin)                                                   1288
      alf=eqs**(1.0/(nlam-1))                                              1288
11481 continue                                                             1289
      rsq=0.0                                                              1289
      a=0.0                                                                1289
      mm=0                                                                 1289
      nlp=0                                                                1289
      nin=nlp                                                              1289
      iz=0                                                                 1289
      mnl=min(mnlam,nlam)                                                  1290
11490 do 11491 m=1,nlam                                                    1291
      if(flmin .lt. 1.0)goto 11511                                         1291
      alm=ulam(m)                                                          1291
      goto 11501                                                           1292
11511 if(m .le. 2)goto 11521                                               1292
      alm=alm*alf                                                          1292
      goto 11501                                                           1293
11521 if(m .ne. 1)goto 11531                                               1293
      alm=big                                                              1293
      goto 11541                                                           1294
11531 continue                                                             1294
      alm=0.0                                                              1295
11550 do 11551 j=1,ni                                                      1295
      if(ju(j).eq.0)goto 11551                                             1295
      if(vp(j).le.0.0)goto 11551                                           1296
      alm=max(alm,abs(g(j))/vp(j))                                         1297
11551 continue                                                             1298
11552 continue                                                             1298
      alm=alf*alm/max(bta,1.0d-3)                                          1299
11541 continue                                                             1300
11501 continue                                                             1300
      dem=alm*omb                                                          1300
      ab=alm*bta                                                           1300
      rsq0=rsq                                                             1300
      jz=1                                                                 1301
11560 continue                                                             1301
11561 continue                                                             1301
      if(iz*jz.ne.0) go to 10360                                           1301
      nlp=nlp+1                                                            1301
      dlx=0.0                                                              1302
11570 do 11571 k=1,ni                                                      1302
      if(ju(k).eq.0)goto 11571                                             1303
      ak=a(k)                                                              1303
      u=g(k)+ak*xv(k)                                                      1303
      v=abs(u)-vp(k)*ab                                                    1303
      a(k)=0.0                                                             1305
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   1306 
     *em)))
      if(a(k).eq.ak)goto 11571                                             1307
      if(mm(k) .ne. 0)goto 11591                                           1307
      nin=nin+1                                                            1307
      if(nin.gt.nx)goto 11572                                              1308
11600 do 11601 j=1,ni                                                      1308
      if(ju(j).eq.0)goto 11601                                             1309
      if(mm(j) .eq. 0)goto 11621                                           1309
      c(j,nin)=c(k,mm(j))                                                  1309
      goto 11601                                                           1309
11621 continue                                                             1310
      if(j .ne. k)goto 11641                                               1310
      c(j,nin)=xv(j)                                                       1310
      goto 11601                                                           1310
11641 continue                                                             1311
      c(j,nin)=  (row_prod(j,k,ix,jx,x,w)-xm(j)*xm(k))/(xs(j)*xs(k))       1313
11601 continue                                                             1314
11602 continue                                                             1314
      mm(k)=nin                                                            1314
      ia(nin)=k                                                            1315
11591 continue                                                             1316
      del=a(k)-ak                                                          1316
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                     1317
      dlx=max(xv(k)*del**2,dlx)                                            1318
11650 do 11651 j=1,ni                                                      1318
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del                              1318
11651 continue                                                             1319
11652 continue                                                             1319
11571 continue                                                             1320
11572 continue                                                             1320
      if(dlx.lt.thr)goto 11562                                             1320
      if(nin.gt.nx)goto 11562                                              1321
      if(nlp .le. maxit)goto 11671                                         1321
      jerr=-m                                                              1321
      return                                                               1321
11671 continue                                                             1322
10360 continue                                                             1322
      iz=1                                                                 1322
      da(1:nin)=a(ia(1:nin))                                               1323
11680 continue                                                             1323
11681 continue                                                             1323
      nlp=nlp+1                                                            1323
      dlx=0.0                                                              1324
11690 do 11691 l=1,nin                                                     1324
      k=ia(l)                                                              1325
      ak=a(k)                                                              1325
      u=g(k)+ak*xv(k)                                                      1325
      v=abs(u)-vp(k)*ab                                                    1325
      a(k)=0.0                                                             1327
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   1328 
     *em)))
      if(a(k).eq.ak)goto 11691                                             1329
      del=a(k)-ak                                                          1329
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                     1330
      dlx=max(xv(k)*del**2,dlx)                                            1331
11700 do 11701 j=1,nin                                                     1331
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                                 1331
11701 continue                                                             1332
11702 continue                                                             1332
11691 continue                                                             1333
11692 continue                                                             1333
      if(dlx.lt.thr)goto 11682                                             1333
      if(nlp .le. maxit)goto 11721                                         1333
      jerr=-m                                                              1333
      return                                                               1333
11721 continue                                                             1334
      goto 11681                                                           1335
11682 continue                                                             1335
      da(1:nin)=a(ia(1:nin))-da(1:nin)                                     1336
11730 do 11731 j=1,ni                                                      1336
      if(mm(j).ne.0)goto 11731                                             1337
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))           1338
11731 continue                                                             1339
11732 continue                                                             1339
      jz=0                                                                 1340
      goto 11561                                                           1341
11562 continue                                                             1341
      if(nin .le. nx)goto 11751                                            1341
      jerr=-10000-m                                                        1341
      goto 11492                                                           1341
11751 continue                                                             1342
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                1342
      kin(m)=nin                                                           1343
      rsqo(m)=rsq                                                          1343
      almo(m)=alm                                                          1343
      lmu=m                                                                1344
      if(m.lt.mnl)goto 11491                                               1344
      if(flmin.ge.1.0)goto 11491                                           1345
      me=0                                                                 1345
11760 do 11761 j=1,nin                                                     1345
      if(ao(j,m).ne.0.0) me=me+1                                           1345
11761 continue                                                             1345
11762 continue                                                             1345
      if(me.gt.ne)goto 11492                                               1346
      if(rsq-rsq0.lt.sml*rsq)goto 11492                                    1346
      if(rsq.gt.rsqmax)goto 11492                                          1347
11491 continue                                                             1348
11492 continue                                                             1348
      deallocate(a,mm,c,da)                                                1349
      return                                                               1350
      end                                                                  1351
      subroutine spelnetn(parm,no,ni,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam,flm   1353 
     *in,ulam,  thr,isd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   1354
      double precision x(*),vp(ni),y(no),w(no),ulam(nlam),cl(2,ni)         1355
      double precision ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)            1356
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1357
      double precision, dimension (:), allocatable :: xm,xs,xv,vlam             
      integer, dimension (:), allocatable :: ju                                 
      allocate(xm(1:ni),stat=jerr)                                         1362
      if(jerr.ne.0) return                                                 1363
      allocate(xs(1:ni),stat=jerr)                                         1364
      if(jerr.ne.0) return                                                 1365
      allocate(ju(1:ni),stat=jerr)                                         1366
      if(jerr.ne.0) return                                                 1367
      allocate(xv(1:ni),stat=jerr)                                         1368
      if(jerr.ne.0) return                                                 1369
      allocate(vlam(1:nlam),stat=jerr)                                     1370
      if(jerr.ne.0) return                                                 1371
      call spchkvars(no,ni,x,ix,ju)                                        1372
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1373
      if(maxval(ju) .gt. 0)goto 11781                                      1373
      jerr=7777                                                            1373
      return                                                               1373
11781 continue                                                             1374
      call spstandard1(no,ni,x,ix,jx,y,w,ju,isd,intr,xm,xs,ym,ys,xv,jerr   1375 
     *)
      if(jerr.ne.0) return                                                 1376
      cl=cl/ys                                                             1376
      if(isd .le. 0)goto 11801                                             1376
11810 do 11811 j=1,ni                                                      1376
      cl(:,j)=cl(:,j)*xs(j)                                                1376
11811 continue                                                             1376
11812 continue                                                             1376
11801 continue                                                             1377
      if(flmin.ge.1.0) vlam=ulam/ys                                        1378
      call spelnet2(parm,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,cl,nlam,flmin,vla   1380 
     *m,thr,maxit,  xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 1381
11820 do 11821 k=1,lmu                                                     1381
      alm(k)=ys*alm(k)                                                     1381
      nk=nin(k)                                                            1382
11830 do 11831 l=1,nk                                                      1382
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                         1382
11831 continue                                                             1382
11832 continue                                                             1382
      a0(k)=0.0                                                            1383
      if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))          1384
11821 continue                                                             1385
11822 continue                                                             1385
      deallocate(xm,xs,ju,xv,vlam)                                         1386
      return                                                               1387
      end                                                                  1388
      subroutine spstandard1 (no,ni,x,ix,jx,y,w,ju,isd,intr,xm,xs,ym,ys,   1389 
     *xv,jerr)
      implicit double precision(a-h,o-z)                                   1390
      double precision x(*),y(no),w(no),xm(ni),xs(ni),xv(ni)               1391
      integer ix(*),jx(*),ju(ni)                                           1392
      w=w/sum(w)                                                           1393
      if(intr .ne. 0)goto 11851                                            1393
      ym=0.0                                                               1394
      ys=sqrt(dot_product(w,y**2)-dot_product(w,y)**2)                     1394
      y=y/ys                                                               1395
11860 do 11861 j=1,ni                                                      1395
      if(ju(j).eq.0)goto 11861                                             1395
      xm(j)=0.0                                                            1395
      jb=ix(j)                                                             1395
      je=ix(j+1)-1                                                         1396
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)                          1397
      if(isd .eq. 0)goto 11881                                             1397
      xbq=dot_product(w(jx(jb:je)),x(jb:je))**2                            1397
      vc=xv(j)-xbq                                                         1398
      xs(j)=sqrt(vc)                                                       1398
      xv(j)=1.0+xbq/vc                                                     1399
      goto 11891                                                           1400
11881 continue                                                             1400
      xs(j)=1.0                                                            1400
11891 continue                                                             1401
11871 continue                                                             1401
11861 continue                                                             1402
11862 continue                                                             1402
      return                                                               1403
11851 continue                                                             1404
11900 do 11901 j=1,ni                                                      1404
      if(ju(j).eq.0)goto 11901                                             1405
      jb=ix(j)                                                             1405
      je=ix(j+1)-1                                                         1405
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1406
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 1407
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       1408
11901 continue                                                             1409
11902 continue                                                             1409
      if(isd .ne. 0)goto 11921                                             1409
      xs=1.0                                                               1409
      goto 11931                                                           1409
11921 continue                                                             1409
      xv=1.0                                                               1409
11931 continue                                                             1410
11911 continue                                                             1410
      ym=dot_product(w,y)                                                  1410
      y=y-ym                                                               1410
      ys=sqrt(dot_product(w,y**2))                                         1410
      y=y/ys                                                               1411
      return                                                               1412
      end                                                                  1413
      subroutine spelnet2(beta,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,cl,nlam,flm   1415 
     *in,ulam,  thr,maxit,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      implicit double precision(a-h,o-z)                                   1416
      double precision y(no),w(no),x(*),vp(ni),ulam(nlam),cl(2,ni)         1417
      double precision ao(nx,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),x   1418 
     *v(ni)
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                          1419
      double precision, dimension (:), allocatable :: a,g                       
      integer, dimension (:), allocatable :: mm,iy                              
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)               1424
      allocate(a(1:ni),stat=jerr)                                          1425
      if(jerr.ne.0) return                                                 1426
      allocate(mm(1:ni),stat=jerr)                                         1427
      if(jerr.ne.0) return                                                 1428
      allocate(g(1:ni),stat=jerr)                                          1429
      if(jerr.ne.0) return                                                 1430
      allocate(iy(1:ni),stat=jerr)                                         1431
      if(jerr.ne.0) return                                                 1432
      bta=beta                                                             1432
      omb=1.0-bta                                                          1432
      alm=0.0                                                              1432
      iy=0                                                                 1434
      alf=1.0                                                              1436
      if(flmin .ge. 1.0)goto 11951                                         1436
      eqs=max(eps,flmin)                                                   1436
      alf=eqs**(1.0/(nlam-1))                                              1436
11951 continue                                                             1437
      rsq=0.0                                                              1437
      a=0.0                                                                1437
      mm=0                                                                 1437
      o=0.0                                                                1437
      nlp=0                                                                1437
      nin=nlp                                                              1437
      iz=0                                                                 1437
      mnl=min(mnlam,nlam)                                                  1438
11960 do 11961 j=1,ni                                                      1438
      if(ju(j).eq.0)goto 11961                                             1439
      jb=ix(j)                                                             1439
      je=ix(j+1)-1                                                         1440
      g(j)=abs(dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(j))    1441
11961 continue                                                             1442
11962 continue                                                             1442
11970 do 11971 m=1,nlam                                                    1442
      alm0=alm                                                             1443
      if(flmin .lt. 1.0)goto 11991                                         1443
      alm=ulam(m)                                                          1443
      goto 11981                                                           1444
11991 if(m .le. 2)goto 12001                                               1444
      alm=alm*alf                                                          1444
      goto 11981                                                           1445
12001 if(m .ne. 1)goto 12011                                               1445
      alm=big                                                              1445
      goto 12021                                                           1446
12011 continue                                                             1446
      alm0=0.0                                                             1447
12030 do 12031 j=1,ni                                                      1447
      if(ju(j).eq.0)goto 12031                                             1447
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           1447
12031 continue                                                             1448
12032 continue                                                             1448
      alm0=alm0/max(bta,1.0d-3)                                            1448
      alm=alf*alm0                                                         1449
12021 continue                                                             1450
11981 continue                                                             1450
      dem=alm*omb                                                          1450
      ab=alm*bta                                                           1450
      rsq0=rsq                                                             1450
      jz=1                                                                 1451
      tlam=bta*(2.0*alm-alm0)                                              1452
12040 do 12041 k=1,ni                                                      1452
      if(iy(k).eq.1)goto 12041                                             1452
      if(ju(k).eq.0)goto 12041                                             1453
      if(g(k).gt.tlam*vp(k)) iy(k)=1                                       1454
12041 continue                                                             1455
12042 continue                                                             1455
12050 continue                                                             1455
12051 continue                                                             1455
      if(iz*jz.ne.0) go to 10360                                           1456
11020 continue                                                             1456
      nlp=nlp+1                                                            1456
      dlx=0.0                                                              1457
12060 do 12061 k=1,ni                                                      1457
      if(iy(k).eq.0)goto 12061                                             1457
      jb=ix(k)                                                             1457
      je=ix(k+1)-1                                                         1458
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)           1459
      ak=a(k)                                                              1459
      u=gk+ak*xv(k)                                                        1459
      v=abs(u)-vp(k)*ab                                                    1459
      a(k)=0.0                                                             1461
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   1462 
     *em)))
      if(a(k).eq.ak)goto 12061                                             1463
      if(mm(k) .ne. 0)goto 12081                                           1463
      nin=nin+1                                                            1463
      if(nin.gt.nx)goto 12062                                              1464
      mm(k)=nin                                                            1464
      ia(nin)=k                                                            1465
12081 continue                                                             1466
      del=a(k)-ak                                                          1466
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1467
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)                         1468
      o=o+del*xm(k)/xs(k)                                                  1468
      dlx=max(xv(k)*del**2,dlx)                                            1469
12061 continue                                                             1470
12062 continue                                                             1470
      if(nin.gt.nx)goto 12052                                              1471
      if(dlx .ge. thr)goto 12101                                           1471
      ixx=0                                                                1472
12110 do 12111 j=1,ni                                                      1472
      if(iy(j).eq.1)goto 12111                                             1472
      if(ju(j).eq.0)goto 12111                                             1473
      jb=ix(j)                                                             1473
      je=ix(j+1)-1                                                         1474
      g(j)=abs(dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(j))    1475
      if(g(j) .le. ab*vp(j))goto 12131                                     1475
      iy(j)=1                                                              1475
      ixx=1                                                                1475
12131 continue                                                             1476
12111 continue                                                             1477
12112 continue                                                             1477
      if(ixx.eq.1) go to 11020                                             1478
      goto 12052                                                           1479
12101 continue                                                             1480
      if(nlp .le. maxit)goto 12151                                         1480
      jerr=-m                                                              1480
      return                                                               1480
12151 continue                                                             1481
10360 continue                                                             1481
      iz=1                                                                 1482
12160 continue                                                             1482
12161 continue                                                             1482
      nlp=nlp+1                                                            1482
      dlx=0.0                                                              1483
12170 do 12171 l=1,nin                                                     1483
      k=ia(l)                                                              1483
      jb=ix(k)                                                             1483
      je=ix(k+1)-1                                                         1484
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)           1485
      ak=a(k)                                                              1485
      u=gk+ak*xv(k)                                                        1485
      v=abs(u)-vp(k)*ab                                                    1485
      a(k)=0.0                                                             1487
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   1488 
     *em)))
      if(a(k).eq.ak)goto 12171                                             1489
      del=a(k)-ak                                                          1489
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1490
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)                         1491
      o=o+del*xm(k)/xs(k)                                                  1491
      dlx=max(xv(k)*del**2,dlx)                                            1492
12171 continue                                                             1493
12172 continue                                                             1493
      if(dlx.lt.thr)goto 12162                                             1493
      if(nlp .le. maxit)goto 12191                                         1493
      jerr=-m                                                              1493
      return                                                               1493
12191 continue                                                             1494
      goto 12161                                                           1495
12162 continue                                                             1495
      jz=0                                                                 1496
      goto 12051                                                           1497
12052 continue                                                             1497
      if(nin .le. nx)goto 12211                                            1497
      jerr=-10000-m                                                        1497
      goto 11972                                                           1497
12211 continue                                                             1498
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                1498
      kin(m)=nin                                                           1499
      rsqo(m)=rsq                                                          1499
      almo(m)=alm                                                          1499
      lmu=m                                                                1500
      if(m.lt.mnl)goto 11971                                               1500
      if(flmin.ge.1.0)goto 11971                                           1501
      me=0                                                                 1501
12220 do 12221 j=1,nin                                                     1501
      if(ao(j,m).ne.0.0) me=me+1                                           1501
12221 continue                                                             1501
12222 continue                                                             1501
      if(me.gt.ne)goto 11972                                               1502
      if(rsq-rsq0.lt.sml*rsq)goto 11972                                    1502
      if(rsq.gt.rsqmax)goto 11972                                          1503
11971 continue                                                             1504
11972 continue                                                             1504
      deallocate(a,mm,g,iy)                                                1505
      return                                                               1506
      end                                                                  1507
      subroutine spchkvars(no,ni,x,ix,ju)                                  1508
      implicit double precision(a-h,o-z)                                   1509
      double precision x(*)                                                1509
      integer ix(*),ju(ni)                                                 1510
12230 do 12231 j=1,ni                                                      1510
      ju(j)=0                                                              1510
      jb=ix(j)                                                             1510
      nj=ix(j+1)-jb                                                        1510
      if(nj.eq.0)goto 12231                                                1511
      je=ix(j+1)-1                                                         1512
      if(nj .ge. no)goto 12251                                             1512
12260 do 12261 i=jb,je                                                     1512
      if(x(i).eq.0.0)goto 12261                                            1512
      ju(j)=1                                                              1512
      goto 12262                                                           1512
12261 continue                                                             1512
12262 continue                                                             1512
      goto 12271                                                           1513
12251 continue                                                             1513
      t=x(jb)                                                              1513
12280 do 12281 i=jb+1,je                                                   1513
      if(x(i).eq.t)goto 12281                                              1513
      ju(j)=1                                                              1513
      goto 12282                                                           1513
12281 continue                                                             1513
12282 continue                                                             1513
12271 continue                                                             1514
12241 continue                                                             1514
12231 continue                                                             1515
12232 continue                                                             1515
      return                                                               1516
      end                                                                  1517
      subroutine cmodval(a0,ca,ia,nin,x,ix,jx,n,f)                         1518
      implicit double precision(a-h,o-z)                                   1519
      double precision ca(*),x(*),f(n)                                     1519
      integer ia(*),ix(*),jx(*)                                            1520
      f=a0                                                                 1521
12290 do 12291 j=1,nin                                                     1521
      k=ia(j)                                                              1521
      kb=ix(k)                                                             1521
      ke=ix(k+1)-1                                                         1522
      f(jx(kb:ke))=f(jx(kb:ke))+ca(j)*x(kb:ke)                             1523
12291 continue                                                             1524
12292 continue                                                             1524
      return                                                               1525
      end                                                                  1526
      function row_prod(i,j,ia,ja,ra,w)                                    1527
      implicit double precision(a-h,o-z)                                   1528
      integer ia(*),ja(*)                                                  1528
      double precision ra(*),w(*)                                          1529
      row_prod=dot(ra(ia(i)),ra(ia(j)),ja(ia(i)),ja(ia(j)),  ia(i+1)-ia(   1531 
     *i),ia(j+1)-ia(j),w)
      return                                                               1532
      end                                                                  1533
      function dot(x,y,mx,my,nx,ny,w)                                      1534
      implicit double precision(a-h,o-z)                                   1535
      double precision x(*),y(*),w(*)                                      1535
      integer mx(*),my(*)                                                  1536
      i=1                                                                  1536
      j=i                                                                  1536
      s=0.0                                                                1537
12300 continue                                                             1537
12301 continue                                                             1537
12310 continue                                                             1538
12311 if(mx(i).ge.my(j))goto 12312                                         1538
      i=i+1                                                                1538
      if(i.gt.nx) go to 12320                                              1538
      goto 12311                                                           1539
12312 continue                                                             1539
      if(mx(i).eq.my(j)) go to 12330                                       1540
12340 continue                                                             1540
12341 if(my(j).ge.mx(i))goto 12342                                         1540
      j=j+1                                                                1540
      if(j.gt.ny) go to 12320                                              1540
      goto 12341                                                           1541
12342 continue                                                             1541
      if(mx(i).eq.my(j)) go to 12330                                       1541
      goto 12301                                                           1542
12330 continue                                                             1542
      s=s+w(mx(i))*x(i)*y(j)                                               1543
      i=i+1                                                                1543
      if(i.gt.nx)goto 12302                                                1543
      j=j+1                                                                1543
      if(j.gt.ny)goto 12302                                                1544
      goto 12301                                                           1545
12302 continue                                                             1545
12320 continue                                                             1545
      dot=s                                                                1546
      return                                                               1547
      end                                                                  1548
      subroutine lognet (parm,no,ni,nc,x,y,g,jd,vp,cl,ne,nx,nlam,flmin,u   1550 
     *lam,thr,  isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,je
     *rr)
      implicit double precision(a-h,o-z)                                   1551
      double precision x(no,ni),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nla   1552 
     *m)
      double precision ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl   1553 
     *(2,ni)
      integer jd(*),ia(nx),nin(nlam)                                       1554
      double precision, dimension (:), allocatable :: xm,xs,ww,vq,xv            
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 12361                                    1558
      jerr=10000                                                           1558
      return                                                               1558
12361 continue                                                             1559
      allocate(ww(1:no),stat=jerr)                                         1560
      if(jerr.ne.0) return                                                 1561
      allocate(ju(1:ni),stat=jerr)                                         1562
      if(jerr.ne.0) return                                                 1563
      allocate(vq(1:ni),stat=jerr)                                         1564
      if(jerr.ne.0) return                                                 1565
      allocate(xm(1:ni),stat=jerr)                                         1566
      if(jerr.ne.0) return                                                 1567
      if(kopt .ne. 2)goto 12381                                            1567
      allocate(xv(1:ni),stat=jerr)                                         1567
      if(jerr.ne.0) return                                                 1567
12381 continue                                                             1568
      if(isd .le. 0)goto 12401                                             1568
      allocate(xs(1:ni),stat=jerr)                                         1568
      if(jerr.ne.0) return                                                 1568
12401 continue                                                             1570
      call chkvars(no,ni,x,ju)                                             1571
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1572
      if(maxval(ju) .gt. 0)goto 12421                                      1572
      jerr=7777                                                            1572
      return                                                               1572
12421 continue                                                             1573
      vq=max(0d0,vp)                                                       1573
      vq=vq*ni/sum(vq)                                                     1574
12430 do 12431 i=1,no                                                      1574
      ww(i)=sum(y(i,:))                                                    1574
      if(ww(i).gt.0.0) y(i,:)=y(i,:)/ww(i)                                 1574
12431 continue                                                             1575
12432 continue                                                             1575
      sw=sum(ww)                                                           1575
      ww=ww/sw                                                             1576
      if(nc .ne. 1)goto 12451                                              1576
      call lstandard1(no,ni,x,ww,ju,isd,intr,xm,xs)                        1577
      if(isd .le. 0)goto 12471                                             1577
12480 do 12481 j=1,ni                                                      1577
      cl(:,j)=cl(:,j)*xs(j)                                                1577
12481 continue                                                             1577
12482 continue                                                             1577
12471 continue                                                             1578
      call lognet2n(parm,no,ni,x,y(:,1),g(:,1),ww,ju,vq,cl,ne,nx,nlam,fl   1580 
     *min,ulam,  thr,isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,n
     *lp,jerr)
      goto 12441                                                           1581
12451 if(kopt .ne. 2)goto 12491                                            1581
      call multlstandard1(no,ni,x,ww,ju,isd,intr,xm,xs,xv)                 1582
      if(isd .le. 0)goto 12511                                             1582
12520 do 12521 j=1,ni                                                      1582
      cl(:,j)=cl(:,j)*xs(j)                                                1582
12521 continue                                                             1582
12522 continue                                                             1582
12511 continue                                                             1583
      call multlognetn(parm,no,ni,nc,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,   1585 
     *ulam,thr,  intr,maxit,xv,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      goto 12531                                                           1586
12491 continue                                                             1586
      call lstandard1(no,ni,x,ww,ju,isd,intr,xm,xs)                        1587
      if(isd .le. 0)goto 12551                                             1587
12560 do 12561 j=1,ni                                                      1587
      cl(:,j)=cl(:,j)*xs(j)                                                1587
12561 continue                                                             1587
12562 continue                                                             1587
12551 continue                                                             1588
      call lognetn(parm,no,ni,nc,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam   1590 
     *,thr,  isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
12531 continue                                                             1591
12441 continue                                                             1591
      if(jerr.gt.0) return                                                 1591
      dev0=2.0*sw*dev0                                                     1592
12570 do 12571 k=1,lmu                                                     1592
      nk=nin(k)                                                            1593
12580 do 12581 ic=1,nc                                                     1593
      if(isd .le. 0)goto 12601                                             1593
12610 do 12611 l=1,nk                                                      1593
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                      1593
12611 continue                                                             1593
12612 continue                                                             1593
12601 continue                                                             1594
      if(intr .ne. 0)goto 12631                                            1594
      a0(ic,k)=0.0                                                         1594
      goto 12641                                                           1595
12631 continue                                                             1595
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))            1595
12641 continue                                                             1596
12621 continue                                                             1596
12581 continue                                                             1597
12582 continue                                                             1597
12571 continue                                                             1598
12572 continue                                                             1598
      deallocate(ww,ju,vq,xm)                                              1598
      if(isd.gt.0) deallocate(xs)                                          1599
      if(kopt.eq.2) deallocate(xv)                                         1600
      return                                                               1601
      end                                                                  1602
      subroutine lstandard1 (no,ni,x,w,ju,isd,intr,xm,xs)                  1603
      implicit double precision(a-h,o-z)                                   1604
      double precision x(no,ni),w(no),xm(ni),xs(ni)                        1604
      integer ju(ni)                                                       1605
      if(intr .ne. 0)goto 12661                                            1606
12670 do 12671 j=1,ni                                                      1606
      if(ju(j).eq.0)goto 12671                                             1606
      xm(j)=0.0                                                            1607
      if(isd .eq. 0)goto 12691                                             1607
      vc=dot_product(w,x(:,j)**2)-dot_product(w,x(:,j))**2                 1608
      xs(j)=sqrt(vc)                                                       1608
      x(:,j)=x(:,j)/xs(j)                                                  1609
12691 continue                                                             1610
12671 continue                                                             1611
12672 continue                                                             1611
      return                                                               1612
12661 continue                                                             1613
12700 do 12701 j=1,ni                                                      1613
      if(ju(j).eq.0)goto 12701                                             1614
      xm(j)=dot_product(w,x(:,j))                                          1614
      x(:,j)=x(:,j)-xm(j)                                                  1615
      if(isd .le. 0)goto 12721                                             1615
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                                 1615
      x(:,j)=x(:,j)/xs(j)                                                  1615
12721 continue                                                             1616
12701 continue                                                             1617
12702 continue                                                             1617
      return                                                               1618
      end                                                                  1619
      subroutine multlstandard1 (no,ni,x,w,ju,isd,intr,xm,xs,xv)           1620
      implicit double precision(a-h,o-z)                                   1621
      double precision x(no,ni),w(no),xm(ni),xs(ni),xv(ni)                 1621
      integer ju(ni)                                                       1622
      if(intr .ne. 0)goto 12741                                            1623
12750 do 12751 j=1,ni                                                      1623
      if(ju(j).eq.0)goto 12751                                             1623
      xm(j)=0.0                                                            1624
      xv(j)=dot_product(w,x(:,j)**2)                                       1625
      if(isd .eq. 0)goto 12771                                             1625
      xbq=dot_product(w,x(:,j))**2                                         1625
      vc=xv(j)-xbq                                                         1626
      xs(j)=sqrt(vc)                                                       1626
      x(:,j)=x(:,j)/xs(j)                                                  1626
      xv(j)=1.0+xbq/vc                                                     1627
12771 continue                                                             1628
12751 continue                                                             1629
12752 continue                                                             1629
      return                                                               1630
12741 continue                                                             1631
12780 do 12781 j=1,ni                                                      1631
      if(ju(j).eq.0)goto 12781                                             1632
      xm(j)=dot_product(w,x(:,j))                                          1632
      x(:,j)=x(:,j)-xm(j)                                                  1633
      xv(j)=dot_product(w,x(:,j)**2)                                       1634
      if(isd .le. 0)goto 12801                                             1634
      xs(j)=sqrt(xv(j))                                                    1634
      x(:,j)=x(:,j)/xs(j)                                                  1634
      xv(j)=1.0                                                            1634
12801 continue                                                             1635
12781 continue                                                             1636
12782 continue                                                             1636
      return                                                               1637
      end                                                                  1638
      subroutine lognet2n(parm,no,ni,x,y,g,w,ju,vp,cl,ne,nx,nlam,flmin,u   1640 
     *lam,shri,  isd,intr,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jer
     *r)
      implicit double precision(a-h,o-z)                                   1641
      double precision x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam),cl(2   1642 
     *,ni)
      double precision a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)             1643
      integer ju(ni),m(nx),kin(nlam)                                       1644
      double precision, dimension (:), allocatable :: b,bs,v,r,xv,q,ga          
      integer, dimension (:), allocatable :: mm,ixx                             
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               1649
      allocate(b(0:ni),stat=jerr)                                          1650
      if(jerr.ne.0) return                                                 1651
      allocate(xv(1:ni),stat=jerr)                                         1652
      if(jerr.ne.0) return                                                 1653
      allocate(ga(1:ni),stat=jerr)                                         1654
      if(jerr.ne.0) return                                                 1655
      allocate(bs(0:ni),stat=jerr)                                         1656
      if(jerr.ne.0) return                                                 1657
      allocate(mm(1:ni),stat=jerr)                                         1658
      if(jerr.ne.0) return                                                 1659
      allocate(ixx(1:ni),stat=jerr)                                        1660
      if(jerr.ne.0) return                                                 1661
      allocate(r(1:no),stat=jerr)                                          1662
      if(jerr.ne.0) return                                                 1663
      allocate(v(1:no),stat=jerr)                                          1664
      if(jerr.ne.0) return                                                 1665
      allocate(q(1:no),stat=jerr)                                          1666
      if(jerr.ne.0) return                                                 1667
      fmax=log(1.0/pmin-1.0)                                               1667
      fmin=-fmax                                                           1667
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                      1668
      bta=parm                                                             1668
      omb=1.0-bta                                                          1669
      q0=dot_product(w,y)                                                  1669
      if(q0 .gt. pmin)goto 12821                                           1669
      jerr=8001                                                            1669
      return                                                               1669
12821 continue                                                             1670
      if(q0 .lt. 1.0-pmin)goto 12841                                       1670
      jerr=9001                                                            1670
      return                                                               1670
12841 continue                                                             1671
      if(intr.eq.0.0) q0=0.5                                               1672
      ixx=0                                                                1672
      al=0.0                                                               1672
      bz=0.0                                                               1672
      if(intr.ne.0) bz=log(q0/(1.0-q0))                                    1673
      if(nonzero(no,g) .ne. 0)goto 12861                                   1673
      vi=q0*(1.0-q0)                                                       1673
      b(0)=bz                                                              1673
      v=vi*w                                                               1674
      r=w*(y-q0)                                                           1674
      q=q0                                                                 1674
      xmz=vi                                                               1674
      dev1=-(bz*q0+log(1.0-q0))                                            1675
      goto 12871                                                           1676
12861 continue                                                             1676
      b(0)=0.0                                                             1677
      if(intr .eq. 0)goto 12891                                            1677
      b(0)=azero(no,y,g,w,jerr)                                            1677
      if(jerr.ne.0) return                                                 1677
12891 continue                                                             1678
      q=1.0/(1.0+exp(-b(0)-g))                                             1678
      v=w*q*(1.0-q)                                                        1678
      r=w*(y-q)                                                            1678
      xmz=sum(v)                                                           1679
      dev1=-(b(0)*q0+dot_product(w,y*g+log(1.0-q)))                        1680
12871 continue                                                             1681
12851 continue                                                             1681
      if(kopt .le. 0)goto 12911                                            1682
      if(isd .le. 0 .or. intr .eq. 0)goto 12931                            1682
      xv=0.25                                                              1682
      goto 12941                                                           1683
12931 continue                                                             1683
12950 do 12951 j=1,ni                                                      1683
      if(ju(j).ne.0) xv(j)=0.25*dot_product(w,x(:,j)**2)                   1683
12951 continue                                                             1683
12952 continue                                                             1683
12941 continue                                                             1684
12921 continue                                                             1684
12911 continue                                                             1685
      dev0=dev1                                                            1686
12960 do 12961 i=1,no                                                      1686
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                        1687
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))              1688
12961 continue                                                             1690
12962 continue                                                             1690
      alf=1.0                                                              1692
      if(flmin .ge. 1.0)goto 12981                                         1692
      eqs=max(eps,flmin)                                                   1692
      alf=eqs**(1.0/(nlam-1))                                              1692
12981 continue                                                             1693
      m=0                                                                  1693
      mm=0                                                                 1693
      nlp=0                                                                1693
      nin=nlp                                                              1693
      mnl=min(mnlam,nlam)                                                  1693
      bs=0.0                                                               1693
      b(1:ni)=0.0                                                          1694
      shr=shri*dev0                                                        1695
12990 do 12991 j=1,ni                                                      1695
      if(ju(j).eq.0)goto 12991                                             1695
      ga(j)=abs(dot_product(r,x(:,j)))                                     1695
12991 continue                                                             1696
12992 continue                                                             1696
13000 do 13001 ilm=1,nlam                                                  1696
      al0=al                                                               1697
      if(flmin .lt. 1.0)goto 13021                                         1697
      al=ulam(ilm)                                                         1697
      goto 13011                                                           1698
13021 if(ilm .le. 2)goto 13031                                             1698
      al=al*alf                                                            1698
      goto 13011                                                           1699
13031 if(ilm .ne. 1)goto 13041                                             1699
      al=big                                                               1699
      goto 13051                                                           1700
13041 continue                                                             1700
      al0=0.0                                                              1701
13060 do 13061 j=1,ni                                                      1701
      if(ju(j).eq.0)goto 13061                                             1701
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            1701
13061 continue                                                             1702
13062 continue                                                             1702
      al0=al0/max(bta,1.0d-3)                                              1702
      al=alf*al0                                                           1703
13051 continue                                                             1704
13011 continue                                                             1704
      al2=al*omb                                                           1704
      al1=al*bta                                                           1704
      tlam=bta*(2.0*al-al0)                                                1705
13070 do 13071 k=1,ni                                                      1705
      if(ixx(k).eq.1)goto 13071                                            1705
      if(ju(k).eq.0)goto 13071                                             1706
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     1707
13071 continue                                                             1708
13072 continue                                                             1708
11020 continue                                                             1709
13080 continue                                                             1709
13081 continue                                                             1709
      bs(0)=b(0)                                                           1709
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                1710
      if(kopt .ne. 0)goto 13101                                            1711
13110 do 13111 j=1,ni                                                      1711
      if(ixx(j).gt.0) xv(j)=dot_product(v,x(:,j)**2)                       1711
13111 continue                                                             1712
13112 continue                                                             1712
13101 continue                                                             1713
13120 continue                                                             1713
13121 continue                                                             1713
      nlp=nlp+1                                                            1713
      dlx=0.0                                                              1714
13130 do 13131 k=1,ni                                                      1714
      if(ixx(k).eq.0)goto 13131                                            1715
      bk=b(k)                                                              1715
      gk=dot_product(r,x(:,k))                                             1716
      u=gk+xv(k)*b(k)                                                      1716
      au=abs(u)-vp(k)*al1                                                  1717
      if(au .gt. 0.0)goto 13151                                            1717
      b(k)=0.0                                                             1717
      goto 13161                                                           1718
13151 continue                                                             1719
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))          1720
13161 continue                                                             1721
13141 continue                                                             1721
      d=b(k)-bk                                                            1721
      if(abs(d).le.0.0)goto 13131                                          1721
      dlx=max(dlx,xv(k)*d**2)                                              1722
      r=r-d*v*x(:,k)                                                       1723
      if(mm(k) .ne. 0)goto 13181                                           1723
      nin=nin+1                                                            1723
      if(nin.gt.nx)goto 13132                                              1724
      mm(k)=nin                                                            1724
      m(nin)=k                                                             1725
13181 continue                                                             1726
13131 continue                                                             1727
13132 continue                                                             1727
      if(nin.gt.nx)goto 13122                                              1728
      d=0.0                                                                1728
      if(intr.ne.0) d=sum(r)/xmz                                           1729
      if(d .eq. 0.0)goto 13201                                             1729
      b(0)=b(0)+d                                                          1729
      dlx=max(dlx,xmz*d**2)                                                1729
      r=r-d*v                                                              1729
13201 continue                                                             1730
      if(dlx.lt.shr)goto 13122                                             1730
      if(nlp .le. maxit)goto 13221                                         1730
      jerr=-ilm                                                            1730
      return                                                               1730
13221 continue                                                             1731
13230 continue                                                             1731
13231 continue                                                             1731
      nlp=nlp+1                                                            1731
      dlx=0.0                                                              1732
13240 do 13241 l=1,nin                                                     1732
      k=m(l)                                                               1732
      bk=b(k)                                                              1733
      gk=dot_product(r,x(:,k))                                             1734
      u=gk+xv(k)*b(k)                                                      1734
      au=abs(u)-vp(k)*al1                                                  1735
      if(au .gt. 0.0)goto 13261                                            1735
      b(k)=0.0                                                             1735
      goto 13271                                                           1736
13261 continue                                                             1737
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))          1738
13271 continue                                                             1739
13251 continue                                                             1739
      d=b(k)-bk                                                            1739
      if(abs(d).le.0.0)goto 13241                                          1739
      dlx=max(dlx,xv(k)*d**2)                                              1740
      r=r-d*v*x(:,k)                                                       1741
13241 continue                                                             1742
13242 continue                                                             1742
      d=0.0                                                                1742
      if(intr.ne.0) d=sum(r)/xmz                                           1743
      if(d .eq. 0.0)goto 13291                                             1743
      b(0)=b(0)+d                                                          1743
      dlx=max(dlx,xmz*d**2)                                                1743
      r=r-d*v                                                              1743
13291 continue                                                             1744
      if(dlx.lt.shr)goto 13232                                             1744
      if(nlp .le. maxit)goto 13311                                         1744
      jerr=-ilm                                                            1744
      return                                                               1744
13311 continue                                                             1745
      goto 13231                                                           1746
13232 continue                                                             1746
      goto 13121                                                           1747
13122 continue                                                             1747
      if(nin.gt.nx)goto 13082                                              1748
13320 do 13321 i=1,no                                                      1748
      fi=b(0)+g(i)                                                         1749
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin)),x(i,m(1:nin)))            1750
      if(fi .ge. fmin)goto 13341                                           1750
      q(i)=0.0                                                             1750
      goto 13331                                                           1750
13341 if(fi .le. fmax)goto 13351                                           1750
      q(i)=1.0                                                             1750
      goto 13361                                                           1751
13351 continue                                                             1751
      q(i)=1.0/(1.0+exp(-fi))                                              1751
13361 continue                                                             1752
13331 continue                                                             1752
13321 continue                                                             1753
13322 continue                                                             1753
      v=w*q*(1.0-q)                                                        1753
      xmz=sum(v)                                                           1753
      if(xmz.le.vmin)goto 13082                                            1753
      r=w*(y-q)                                                            1754
      if(xmz*(b(0)-bs(0))**2 .ge. shr)goto 13381                           1754
      ix=0                                                                 1755
13390 do 13391 j=1,nin                                                     1755
      k=m(j)                                                               1756
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 13391                           1756
      ix=1                                                                 1756
      goto 13392                                                           1757
13391 continue                                                             1758
13392 continue                                                             1758
      if(ix .ne. 0)goto 13411                                              1759
13420 do 13421 k=1,ni                                                      1759
      if(ixx(k).eq.1)goto 13421                                            1759
      if(ju(k).eq.0)goto 13421                                             1760
      ga(k)=abs(dot_product(r,x(:,k)))                                     1761
      if(ga(k) .le. al1*vp(k))goto 13441                                   1761
      ixx(k)=1                                                             1761
      ix=1                                                                 1761
13441 continue                                                             1762
13421 continue                                                             1763
13422 continue                                                             1763
      if(ix.eq.1) go to 11020                                              1764
      goto 13082                                                           1765
13411 continue                                                             1766
13381 continue                                                             1767
      goto 13081                                                           1768
13082 continue                                                             1768
      if(nin .le. nx)goto 13461                                            1768
      jerr=-10000-ilm                                                      1768
      goto 13002                                                           1768
13461 continue                                                             1769
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                1769
      kin(ilm)=nin                                                         1770
      a0(ilm)=b(0)                                                         1770
      alm(ilm)=al                                                          1770
      lmu=ilm                                                              1771
      devi=dev2(no,w,y,q,pmin)                                             1772
      dev(ilm)=(dev1-devi)/dev0                                            1772
      if(xmz.le.vmin)goto 13002                                            1773
      if(ilm.lt.mnl)goto 13001                                             1773
      if(flmin.ge.1.0)goto 13001                                           1774
      me=0                                                                 1774
13470 do 13471 j=1,nin                                                     1774
      if(a(j,ilm).ne.0.0) me=me+1                                          1774
13471 continue                                                             1774
13472 continue                                                             1774
      if(me.gt.ne)goto 13002                                               1775
      if(dev(ilm).gt.devmax)goto 13002                                     1775
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 13002                             1776
13001 continue                                                             1777
13002 continue                                                             1777
      g=log(q/(1.0-q))                                                     1778
      deallocate(b,bs,v,r,xv,q,mm,ga,ixx)                                  1779
      return                                                               1780
      end                                                                  1781
      function dev2(n,w,y,p,pmin)                                          1782
      implicit double precision(a-h,o-z)                                   1783
      double precision w(n),y(n),p(n)                                      1784
      pmax=1.0-pmin                                                        1784
      s=0.0                                                                1785
13480 do 13481 i=1,n                                                       1785
      pi=min(max(pmin,p(i)),pmax)                                          1786
      s=s-w(i)*(y(i)*log(pi)+(1.0-y(i))*log(1.0-pi))                       1787
13481 continue                                                             1788
13482 continue                                                             1788
      dev2=s                                                               1789
      return                                                               1790
      end                                                                  1791
      function azero(n,y,g,q,jerr)                                         1792
      implicit double precision(a-h,o-z)                                   1793
      parameter(eps=1.0d-7)                                                1794
      double precision y(n),g(n),q(n)                                      1795
      double precision, dimension (:), allocatable :: e,p,w                     
      azero = 0.0                                                          1799
      allocate(e(1:n),stat=jerr)                                           1800
      if(jerr.ne.0) return                                                 1801
      allocate(p(1:n),stat=jerr)                                           1802
      if(jerr.ne.0) return                                                 1803
      allocate(w(1:n),stat=jerr)                                           1804
      if(jerr.ne.0) return                                                 1805
      az=0.0                                                               1805
      e=exp(-g)                                                            1805
      qy=dot_product(q,y)                                                  1805
      p=1.0/(1.0+e)                                                        1806
13490 continue                                                             1806
13491 continue                                                             1806
      w=q*p*(1.0-p)                                                        1807
      d=(qy-dot_product(q,p))/sum(w)                                       1807
      az=az+d                                                              1807
      if(abs(d).lt.eps)goto 13492                                          1808
      ea0=exp(-az)                                                         1808
      p=1.0/(1.0+ea0*e)                                                    1809
      goto 13491                                                           1810
13492 continue                                                             1810
      azero=az                                                             1811
      deallocate(e,p,w)                                                    1812
      return                                                               1813
      end                                                                  1814
      subroutine lognetn(parm,no,ni,nc,x,y,g,w,ju,vp,cl,ne,nx,nlam,flmin   1816 
     *,ulam,shri,  isd,intr,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,j
     *err)
      implicit double precision(a-h,o-z)                                   1817
      double precision x(no,ni),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam   1818 
     *)
      double precision a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(   1819 
     *2,ni)
      integer ju(ni),m(nx),kin(nlam)                                       1820
      double precision, dimension (:,:), allocatable :: q                       
      double precision, dimension (:), allocatable :: sxp,sxpl                  
      double precision, dimension (:), allocatable :: di,v,r,ga                 
      double precision, dimension (:,:), allocatable :: b,bs,xv                 
      integer, dimension (:), allocatable :: mm,is,ixx                          
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      if(jerr.ne.0) return                                                      
      allocate(xv(1:ni,1:nc),stat=jerr)                                         
      if(jerr.ne.0) return					                                                 
      allocate(bs(0:ni,1:nc),stat=jerr)                                         
      if(jerr.ne.0) return					                                                 
      allocate(q(1:no,1:nc),stat=jerr)                                          
      if(jerr.ne.0) return		                                                    
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               1835
      exmn=-exmx                                                           1836
      allocate(r(1:no),stat=jerr)                                          1837
      if(jerr.ne.0) return                                                 1838
      allocate(v(1:no),stat=jerr)                                          1839
      if(jerr.ne.0) return                                                 1840
      allocate(mm(1:ni),stat=jerr)                                         1841
      if(jerr.ne.0) return                                                 1842
      allocate(is(1:max(nc,ni)),stat=jerr)                                 1843
      if(jerr.ne.0) return                                                 1844
      allocate(sxp(1:no),stat=jerr)                                        1845
      if(jerr.ne.0) return                                                 1846
      allocate(sxpl(1:no),stat=jerr)                                       1847
      if(jerr.ne.0) return                                                 1848
      allocate(di(1:no),stat=jerr)                                         1849
      if(jerr.ne.0) return                                                 1850
      allocate(ga(1:ni),stat=jerr)                                         1851
      if(jerr.ne.0) return                                                 1852
      allocate(ixx(1:ni),stat=jerr)                                        1853
      if(jerr.ne.0) return                                                 1854
      pmax=1.0-pmin                                                        1854
      emin=pmin/pmax                                                       1854
      emax=1.0/emin                                                        1855
      pfm=(1.0+pmin)*pmin                                                  1855
      pfx=(1.0-pmin)*pmax                                                  1855
      vmin=pfm*pmax                                                        1856
      bta=parm                                                             1856
      omb=1.0-bta                                                          1856
      dev1=0.0                                                             1856
      dev0=0.0                                                             1857
13500 do 13501 ic=1,nc                                                     1857
      q0=dot_product(w,y(:,ic))                                            1858
      if(q0 .gt. pmin)goto 13521                                           1858
      jerr =8000+ic                                                        1858
      return                                                               1858
13521 continue                                                             1859
      if(q0 .lt. 1.0-pmin)goto 13541                                       1859
      jerr =9000+ic                                                        1859
      return                                                               1859
13541 continue                                                             1860
      if(intr .ne. 0)goto 13561                                            1860
      q0=1.0/nc                                                            1860
      b(0,ic)=0.0                                                          1860
      goto 13571                                                           1861
13561 continue                                                             1861
      b(0,ic)=log(q0)                                                      1861
      dev1=dev1-q0*b(0,ic)                                                 1861
13571 continue                                                             1862
13551 continue                                                             1862
      b(1:ni,ic)=0.0                                                       1863
13501 continue                                                             1864
13502 continue                                                             1864
      if(intr.eq.0) dev1=log(float(nc))                                    1864
      ixx=0                                                                1864
      al=0.0                                                               1865
      if(nonzero(no*nc,g) .ne. 0)goto 13591                                1866
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         1866
      sxp=0.0                                                              1867
13600 do 13601 ic=1,nc                                                     1867
      q(:,ic)=exp(b(0,ic))                                                 1867
      sxp=sxp+q(:,ic)                                                      1867
13601 continue                                                             1868
13602 continue                                                             1868
      goto 13611                                                           1869
13591 continue                                                             1869
13620 do 13621 i=1,no                                                      1869
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1869
13621 continue                                                             1869
13622 continue                                                             1869
      sxp=0.0                                                              1870
      if(intr .ne. 0)goto 13641                                            1870
      b(0,:)=0.0                                                           1870
      goto 13651                                                           1871
13641 continue                                                             1871
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 1871
      if(jerr.ne.0) return                                                 1871
13651 continue                                                             1872
13631 continue                                                             1872
      dev1=0.0                                                             1873
13660 do 13661 ic=1,nc                                                     1873
      q(:,ic)=b(0,ic)+g(:,ic)                                              1874
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             1875
      q(:,ic)=exp(q(:,ic))                                                 1875
      sxp=sxp+q(:,ic)                                                      1876
13661 continue                                                             1877
13662 continue                                                             1877
      sxpl=w*log(sxp)                                                      1877
13670 do 13671 ic=1,nc                                                     1877
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  1877
13671 continue                                                             1878
13672 continue                                                             1878
13611 continue                                                             1879
13581 continue                                                             1879
13680 do 13681 ic=1,nc                                                     1879
13690 do 13691 i=1,no                                                      1879
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               1879
13691 continue                                                             1879
13692 continue                                                             1879
13681 continue                                                             1880
13682 continue                                                             1880
      dev0=dev0+dev1                                                       1881
      if(kopt .le. 0)goto 13711                                            1882
      if(isd .le. 0 .or. intr .eq. 0)goto 13731                            1882
      xv=0.25                                                              1882
      goto 13741                                                           1883
13731 continue                                                             1883
13750 do 13751 j=1,ni                                                      1883
      if(ju(j).ne.0) xv(j,:)=0.25*dot_product(w,x(:,j)**2)                 1883
13751 continue                                                             1883
13752 continue                                                             1883
13741 continue                                                             1884
13721 continue                                                             1884
13711 continue                                                             1886
      alf=1.0                                                              1888
      if(flmin .ge. 1.0)goto 13771                                         1888
      eqs=max(eps,flmin)                                                   1888
      alf=eqs**(1.0/(nlam-1))                                              1888
13771 continue                                                             1889
      m=0                                                                  1889
      mm=0                                                                 1889
      nin=0                                                                1889
      nlp=0                                                                1889
      mnl=min(mnlam,nlam)                                                  1889
      bs=0.0                                                               1889
      shr=shri*dev0                                                        1890
      ga=0.0                                                               1891
13780 do 13781 ic=1,nc                                                     1891
      r=w*(y(:,ic)-q(:,ic)/sxp)                                            1892
13790 do 13791 j=1,ni                                                      1892
      if(ju(j).ne.0) ga(j)=max(ga(j),abs(dot_product(r,x(:,j))))           1892
13791 continue                                                             1893
13792 continue                                                             1893
13781 continue                                                             1894
13782 continue                                                             1894
13800 do 13801 ilm=1,nlam                                                  1894
      al0=al                                                               1895
      if(flmin .lt. 1.0)goto 13821                                         1895
      al=ulam(ilm)                                                         1895
      goto 13811                                                           1896
13821 if(ilm .le. 2)goto 13831                                             1896
      al=al*alf                                                            1896
      goto 13811                                                           1897
13831 if(ilm .ne. 1)goto 13841                                             1897
      al=big                                                               1897
      goto 13851                                                           1898
13841 continue                                                             1898
      al0=0.0                                                              1899
13860 do 13861 j=1,ni                                                      1899
      if(ju(j).eq.0)goto 13861                                             1899
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            1899
13861 continue                                                             1900
13862 continue                                                             1900
      al0=al0/max(bta,1.0d-3)                                              1900
      al=alf*al0                                                           1901
13851 continue                                                             1902
13811 continue                                                             1902
      al2=al*omb                                                           1902
      al1=al*bta                                                           1902
      tlam=bta*(2.0*al-al0)                                                1903
13870 do 13871 k=1,ni                                                      1903
      if(ixx(k).eq.1)goto 13871                                            1903
      if(ju(k).eq.0)goto 13871                                             1904
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     1905
13871 continue                                                             1906
13872 continue                                                             1906
11020 continue                                                             1907
13880 continue                                                             1907
13881 continue                                                             1907
      ix=0                                                                 1907
      jx=ix                                                                1907
      ig=0                                                                 1908
13890 do 13891 ic=1,nc                                                     1908
      bs(0,ic)=b(0,ic)                                                     1909
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          1910
      xmz=0.0                                                              1911
13900 do 13901 i=1,no                                                      1911
      pic=q(i,ic)/sxp(i)                                                   1912
      if(pic .ge. pfm)goto 13921                                           1912
      pic=0.0                                                              1912
      v(i)=0.0                                                             1912
      goto 13911                                                           1913
13921 if(pic .le. pfx)goto 13931                                           1913
      pic=1.0                                                              1913
      v(i)=0.0                                                             1913
      goto 13941                                                           1914
13931 continue                                                             1914
      v(i)=w(i)*pic*(1.0-pic)                                              1914
      xmz=xmz+v(i)                                                         1914
13941 continue                                                             1915
13911 continue                                                             1915
      r(i)=w(i)*(y(i,ic)-pic)                                              1916
13901 continue                                                             1917
13902 continue                                                             1917
      if(xmz.le.vmin)goto 13891                                            1917
      ig=1                                                                 1918
      if(kopt .ne. 0)goto 13961                                            1919
13970 do 13971 j=1,ni                                                      1919
      if(ixx(j).gt.0) xv(j,ic)=dot_product(v,x(:,j)**2)                    1919
13971 continue                                                             1920
13972 continue                                                             1920
13961 continue                                                             1921
13980 continue                                                             1921
13981 continue                                                             1921
      nlp=nlp+1                                                            1921
      dlx=0.0                                                              1922
13990 do 13991 k=1,ni                                                      1922
      if(ixx(k).eq.0)goto 13991                                            1923
      bk=b(k,ic)                                                           1923
      gk=dot_product(r,x(:,k))                                             1924
      u=gk+xv(k,ic)*b(k,ic)                                                1924
      au=abs(u)-vp(k)*al1                                                  1925
      if(au .gt. 0.0)goto 14011                                            1925
      b(k,ic)=0.0                                                          1925
      goto 14021                                                           1926
14011 continue                                                             1927
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))   1929 
     *)
14021 continue                                                             1930
14001 continue                                                             1930
      d=b(k,ic)-bk                                                         1930
      if(abs(d).le.0.0)goto 13991                                          1931
      dlx=max(dlx,xv(k,ic)*d**2)                                           1931
      r=r-d*v*x(:,k)                                                       1932
      if(mm(k) .ne. 0)goto 14041                                           1932
      nin=nin+1                                                            1933
      if(nin .le. nx)goto 14061                                            1933
      jx=1                                                                 1933
      goto 13992                                                           1933
14061 continue                                                             1934
      mm(k)=nin                                                            1934
      m(nin)=k                                                             1935
14041 continue                                                             1936
13991 continue                                                             1937
13992 continue                                                             1937
      if(jx.gt.0)goto 13982                                                1938
      d=0.0                                                                1938
      if(intr.ne.0) d=sum(r)/xmz                                           1939
      if(d .eq. 0.0)goto 14081                                             1939
      b(0,ic)=b(0,ic)+d                                                    1939
      dlx=max(dlx,xmz*d**2)                                                1939
      r=r-d*v                                                              1939
14081 continue                                                             1940
      if(dlx.lt.shr)goto 13982                                             1941
      if(nlp .le. maxit)goto 14101                                         1941
      jerr=-ilm                                                            1941
      return                                                               1941
14101 continue                                                             1942
14110 continue                                                             1942
14111 continue                                                             1942
      nlp=nlp+1                                                            1942
      dlx=0.0                                                              1943
14120 do 14121 l=1,nin                                                     1943
      k=m(l)                                                               1943
      bk=b(k,ic)                                                           1944
      gk=dot_product(r,x(:,k))                                             1945
      u=gk+xv(k,ic)*b(k,ic)                                                1945
      au=abs(u)-vp(k)*al1                                                  1946
      if(au .gt. 0.0)goto 14141                                            1946
      b(k,ic)=0.0                                                          1946
      goto 14151                                                           1947
14141 continue                                                             1948
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))   1950 
     *)
14151 continue                                                             1951
14131 continue                                                             1951
      d=b(k,ic)-bk                                                         1951
      if(abs(d).le.0.0)goto 14121                                          1952
      dlx=max(dlx,xv(k,ic)*d**2)                                           1952
      r=r-d*v*x(:,k)                                                       1953
14121 continue                                                             1954
14122 continue                                                             1954
      d=0.0                                                                1954
      if(intr.ne.0) d=sum(r)/xmz                                           1955
      if(d .eq. 0.0)goto 14171                                             1955
      b(0,ic)=b(0,ic)+d                                                    1956
      dlx=max(dlx,xmz*d**2)                                                1956
      r=r-d*v                                                              1957
14171 continue                                                             1958
      if(dlx.lt.shr)goto 14112                                             1958
      if(nlp .le. maxit)goto 14191                                         1958
      jerr=-ilm                                                            1958
      return                                                               1958
14191 continue                                                             1959
      goto 14111                                                           1960
14112 continue                                                             1960
      goto 13981                                                           1961
13982 continue                                                             1961
      if(jx.gt.0)goto 13892                                                1962
      if(xmz*(b(0,ic)-bs(0,ic))**2.gt.shr) ix=1                            1963
      if(ix .ne. 0)goto 14211                                              1964
14220 do 14221 j=1,nin                                                     1964
      k=m(j)                                                               1965
      if(xv(k,ic)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 14241                1965
      ix=1                                                                 1965
      goto 14222                                                           1965
14241 continue                                                             1966
14221 continue                                                             1967
14222 continue                                                             1967
14211 continue                                                             1968
14250 do 14251 i=1,no                                                      1968
      fi=b(0,ic)+g(i,ic)                                                   1970
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)))         1971
      fi=min(max(exmn,fi),exmx)                                            1971
      sxp(i)=sxp(i)-q(i,ic)                                                1972
      q(i,ic)=min(max(emin*sxp(i),exp(fi)),emax*sxp(i))                    1973
      sxp(i)=sxp(i)+q(i,ic)                                                1974
14251 continue                                                             1975
14252 continue                                                             1975
13891 continue                                                             1976
13892 continue                                                             1976
      s=-sum(b(0,:))/nc                                                    1976
      b(0,:)=b(0,:)+s                                                      1976
      di=s                                                                 1977
14260 do 14261 j=1,nin                                                     1977
      l=m(j)                                                               1978
      if(vp(l) .gt. 0.0)goto 14281                                         1978
      s=sum(b(l,:))/nc                                                     1978
      goto 14291                                                           1979
14281 continue                                                             1979
      s=elc(parm,nc,cl(:,l),b(l,:),is)                                     1979
14291 continue                                                             1980
14271 continue                                                             1980
      b(l,:)=b(l,:)-s                                                      1980
      di=di-s*x(:,l)                                                       1981
14261 continue                                                             1982
14262 continue                                                             1982
      di=exp(di)                                                           1982
      sxp=sxp*di                                                           1982
14300 do 14301 ic=1,nc                                                     1982
      q(:,ic)=q(:,ic)*di                                                   1982
14301 continue                                                             1983
14302 continue                                                             1983
      if(jx.gt.0)goto 13882                                                1983
      if(ig.eq.0)goto 13882                                                1984
      if(ix .ne. 0)goto 14321                                              1985
14330 do 14331 k=1,ni                                                      1985
      if(ixx(k).eq.1)goto 14331                                            1985
      if(ju(k).eq.0)goto 14331                                             1985
      ga(k)=0.0                                                            1985
14331 continue                                                             1986
14332 continue                                                             1986
14340 do 14341 ic=1,nc                                                     1986
      r=w*(y(:,ic)-q(:,ic)/sxp)                                            1987
14350 do 14351 k=1,ni                                                      1987
      if(ixx(k).eq.1)goto 14351                                            1987
      if(ju(k).eq.0)goto 14351                                             1988
      ga(k)=max(ga(k),abs(dot_product(r,x(:,k))))                          1989
14351 continue                                                             1990
14352 continue                                                             1990
14341 continue                                                             1991
14342 continue                                                             1991
14360 do 14361 k=1,ni                                                      1991
      if(ixx(k).eq.1)goto 14361                                            1991
      if(ju(k).eq.0)goto 14361                                             1992
      if(ga(k) .le. al1*vp(k))goto 14381                                   1992
      ixx(k)=1                                                             1992
      ix=1                                                                 1992
14381 continue                                                             1993
14361 continue                                                             1994
14362 continue                                                             1994
      if(ix.eq.1) go to 11020                                              1995
      goto 13882                                                           1996
14321 continue                                                             1997
      goto 13881                                                           1998
13882 continue                                                             1998
      if(jx .le. 0)goto 14401                                              1998
      jerr=-10000-ilm                                                      1998
      goto 13802                                                           1998
14401 continue                                                             1998
      devi=0.0                                                             1999
14410 do 14411 ic=1,nc                                                     2000
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          2000
      a0(ic,ilm)=b(0,ic)                                                   2001
14420 do 14421 i=1,no                                                      2001
      if(y(i,ic).le.0.0)goto 14421                                         2002
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           2003
14421 continue                                                             2004
14422 continue                                                             2004
14411 continue                                                             2005
14412 continue                                                             2005
      kin(ilm)=nin                                                         2005
      alm(ilm)=al                                                          2005
      lmu=ilm                                                              2006
      dev(ilm)=(dev1-devi)/dev0                                            2006
      if(ig.eq.0)goto 13802                                                2007
      if(ilm.lt.mnl)goto 13801                                             2007
      if(flmin.ge.1.0)goto 13801                                           2008
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 13802             2009
      if(dev(ilm).gt.devmax)goto 13802                                     2009
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 13802                             2010
13801 continue                                                             2011
13802 continue                                                             2011
      g=log(q)                                                             2011
14430 do 14431 i=1,no                                                      2011
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         2011
14431 continue                                                             2012
14432 continue                                                             2012
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,ga,ixx)                           2013
      return                                                               2014
      end                                                                  2015
      subroutine kazero(kk,n,y,g,q,az,jerr)                                2016
      implicit double precision(a-h,o-z)                                   2017
      parameter(eps=1.0d-7)                                                2018
      double precision y(n,kk),g(n,kk),q(n),az(kk)                         2019
      double precision, dimension (:), allocatable :: s                         
      double precision, dimension (:,:), allocatable :: e                       
      allocate(e(1:n,1:kk),stat=jerr)                                           
      if(jerr.ne.0) return                                                      
      allocate(s(1:n),stat=jerr)                                           2026
      if(jerr.ne.0) return                                                 2027
      az=0.0                                                               2027
      e=exp(g)                                                             2027
14440 do 14441 i=1,n                                                       2027
      s(i)=sum(e(i,:))                                                     2027
14441 continue                                                             2028
14442 continue                                                             2028
14450 continue                                                             2028
14451 continue                                                             2028
      dm=0.0                                                               2029
14460 do 14461 k=1,kk                                                      2029
      t=0.0                                                                2029
      u=t                                                                  2030
14470 do 14471 i=1,n                                                       2030
      pik=e(i,k)/s(i)                                                      2031
      t=t+q(i)*(y(i,k)-pik)                                                2031
      u=u+q(i)*pik*(1.0-pik)                                               2032
14471 continue                                                             2033
14472 continue                                                             2033
      d=t/u                                                                2033
      az(k)=az(k)+d                                                        2033
      ed=exp(d)                                                            2033
      dm=max(dm,abs(d))                                                    2034
14480 do 14481 i=1,n                                                       2034
      z=e(i,k)                                                             2034
      e(i,k)=z*ed                                                          2034
      s(i)=s(i)-z+e(i,k)                                                   2034
14481 continue                                                             2035
14482 continue                                                             2035
14461 continue                                                             2036
14462 continue                                                             2036
      if(dm.lt.eps)goto 14452                                              2036
      goto 14451                                                           2037
14452 continue                                                             2037
      az=az-sum(az)/kk                                                     2038
      deallocate(e,s)                                                      2039
      return                                                               2040
      end                                                                  2041
      function elc(parm,n,cl,a,m)                                          2042
      implicit double precision(a-h,o-z)                                   2043
      double precision a(n),cl(2)                                          2043
      integer m(n)                                                         2044
      fn=n                                                                 2044
      am=sum(a)/fn                                                         2045
      if((parm .ne. 0.0) .and. (n .ne. 2))goto 14501                       2045
      elc=am                                                               2045
      go to 14510                                                          2045
14501 continue                                                             2046
14520 do 14521 i=1,n                                                       2046
      m(i)=i                                                               2046
14521 continue                                                             2046
14522 continue                                                             2046
      call psort7(a,m,1,n)                                                 2047
      if(a(m(1)) .ne. a(m(n)))goto 14541                                   2047
      elc=a(1)                                                             2047
      go to 14510                                                          2047
14541 continue                                                             2048
      if(mod(n,2) .ne. 1)goto 14561                                        2048
      ad=a(m(n/2+1))                                                       2048
      goto 14571                                                           2049
14561 continue                                                             2049
      ad=0.5*(a(m(n/2+1))+a(m(n/2)))                                       2049
14571 continue                                                             2050
14551 continue                                                             2050
      if(parm .ne. 1.0)goto 14591                                          2050
      elc=ad                                                               2050
      go to 14510                                                          2050
14591 continue                                                             2051
      b1=min(am,ad)                                                        2051
      b2=max(am,ad)                                                        2051
      k2=1                                                                 2052
14600 continue                                                             2052
14601 if(a(m(k2)).gt.b1)goto 14602                                         2052
      k2=k2+1                                                              2052
      goto 14601                                                           2052
14602 continue                                                             2052
      k1=k2-1                                                              2053
14610 continue                                                             2053
14611 if(a(m(k2)).ge.b2)goto 14612                                         2053
      k2=k2+1                                                              2053
      goto 14611                                                           2054
14612 continue                                                             2054
      r=parm/((1.0-parm)*fn)                                               2054
      is=0                                                                 2054
      sm=n-2*(k1-1)                                                        2055
14620 do 14621 k=k1,k2-1                                                   2055
      sm=sm-2.0                                                            2055
      s=r*sm+am                                                            2056
      if(s .le. a(m(k)) .or. s .gt. a(m(k+1)))goto 14641                   2056
      is=k                                                                 2056
      goto 14622                                                           2056
14641 continue                                                             2057
14621 continue                                                             2058
14622 continue                                                             2058
      if(is .eq. 0)goto 14661                                              2058
      elc=s                                                                2058
      go to 14510                                                          2058
14661 continue                                                             2058
      r2=2.0*r                                                             2058
      s1=a(m(k1))                                                          2058
      am2=2.0*am                                                           2059
      cri=r2*sum(abs(a-s1))+s1*(s1-am2)                                    2059
      elc=s1                                                               2060
14670 do 14671 k=k1+1,k2                                                   2060
      s=a(m(k))                                                            2060
      if(s.eq.s1)goto 14671                                                2061
      c=r2*sum(abs(a-s))+s*(s-am2)                                         2062
      if(c .ge. cri)goto 14691                                             2062
      cri=c                                                                2062
      elc=s                                                                2062
14691 continue                                                             2062
      s1=s                                                                 2063
14671 continue                                                             2064
14672 continue                                                             2064
14510 continue                                                             2064
      elc=max(maxval(a-cl(2)),min(minval(a-cl(1)),elc))                    2065
      return                                                               2066
      end                                                                  2067
      function nintot(ni,nx,nc,a,m,nin,is)                                 2068
      implicit double precision(a-h,o-z)                                   2069
      double precision a(nx,nc)                                            2069
      integer m(nx),is(ni)                                                 2070
      is=0                                                                 2070
      nintot=0                                                             2071
14700 do 14701 ic=1,nc                                                     2071
14710 do 14711 j=1,nin                                                     2071
      k=m(j)                                                               2071
      if(is(k).ne.0)goto 14711                                             2072
      if(a(j,ic).eq.0.0)goto 14711                                         2072
      is(k)=k                                                              2072
      nintot=nintot+1                                                      2073
14711 continue                                                             2073
14712 continue                                                             2073
14701 continue                                                             2074
14702 continue                                                             2074
      return                                                               2075
      end                                                                  2076
      subroutine luncomp(ni,nx,nc,ca,ia,nin,a)                             2077
      implicit double precision(a-h,o-z)                                   2078
      double precision ca(nx,nc),a(ni,nc)                                  2078
      integer ia(nx)                                                       2079
      a=0.0                                                                2080
14720 do 14721 ic=1,nc                                                     2080
      if(nin.gt.0) a(ia(1:nin),ic)=ca(1:nin,ic)                            2080
14721 continue                                                             2081
14722 continue                                                             2081
      return                                                               2082
      end                                                                  2083
      subroutine lmodval(nt,x,nc,nx,a0,ca,ia,nin,ans)                      2084
      implicit double precision(a-h,o-z)                                   2085
      double precision a0(nc),ca(nx,nc),x(nt,*),ans(nc,nt)                 2085
      integer ia(nx)                                                       2086
14730 do 14731 i=1,nt                                                      2086
14740 do 14741 ic=1,nc                                                     2086
      ans(ic,i)=a0(ic)                                                     2088
      if(nin.gt.0) ans(ic,i)=ans(ic,i)+dot_product(ca(1:nin,ic),x(i,ia(1   2089 
     *:nin)))
14741 continue                                                             2089
14742 continue                                                             2089
14731 continue                                                             2090
14732 continue                                                             2090
      return                                                               2091
      end                                                                  2092
      subroutine splognet (parm,no,ni,nc,x,ix,jx,y,g,jd,vp,cl,ne,nx,nlam   2094 
     *,flmin,  ulam,thr,isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,al
     *m,nlp,jerr)
      implicit double precision(a-h,o-z)                                   2095
      double precision x(*),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nlam)     2096
      double precision ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl   2097 
     *(2,ni)
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           2098
      double precision, dimension (:), allocatable :: xm,xs,ww,vq,xv            
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 14761                                    2102
      jerr=10000                                                           2102
      return                                                               2102
14761 continue                                                             2103
      allocate(ww(1:no),stat=jerr)                                         2104
      if(jerr.ne.0) return                                                 2105
      allocate(ju(1:ni),stat=jerr)                                         2106
      if(jerr.ne.0) return                                                 2107
      allocate(vq(1:ni),stat=jerr)                                         2108
      if(jerr.ne.0) return                                                 2109
      allocate(xm(1:ni),stat=jerr)                                         2110
      if(jerr.ne.0) return                                                 2111
      allocate(xs(1:ni),stat=jerr)                                         2112
      if(jerr.ne.0) return                                                 2113
      if(kopt .ne. 2)goto 14781                                            2113
      allocate(xv(1:ni),stat=jerr)                                         2113
      if(jerr.ne.0) return                                                 2113
14781 continue                                                             2115
      call spchkvars(no,ni,x,ix,ju)                                        2116
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2117
      if(maxval(ju) .gt. 0)goto 14801                                      2117
      jerr=7777                                                            2117
      return                                                               2117
14801 continue                                                             2118
      vq=max(0d0,vp)                                                       2118
      vq=vq*ni/sum(vq)                                                     2119
14810 do 14811 i=1,no                                                      2119
      ww(i)=sum(y(i,:))                                                    2119
      if(ww(i).gt.0.0) y(i,:)=y(i,:)/ww(i)                                 2119
14811 continue                                                             2120
14812 continue                                                             2120
      sw=sum(ww)                                                           2120
      ww=ww/sw                                                             2121
      if(nc .ne. 1)goto 14831                                              2121
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs)                2122
      if(isd .le. 0)goto 14851                                             2122
14860 do 14861 j=1,ni                                                      2122
      cl(:,j)=cl(:,j)*xs(j)                                                2122
14861 continue                                                             2122
14862 continue                                                             2122
14851 continue                                                             2123
      call sprlognet2n(parm,no,ni,x,ix,jx,y(:,1),g(:,1),ww,ju,vq,cl,ne,n   2126 
     *x,nlam,  flmin,ulam,thr,isd,intr,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin
     *,dev0,dev,  alm,nlp,jerr)
      goto 14821                                                           2127
14831 if(kopt .ne. 2)goto 14871                                            2128
      call multsplstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs,xv)         2129
      if(isd .le. 0)goto 14891                                             2129
14900 do 14901 j=1,ni                                                      2129
      cl(:,j)=cl(:,j)*xs(j)                                                2129
14901 continue                                                             2129
14902 continue                                                             2129
14891 continue                                                             2130
      call multsprlognetn(parm,no,ni,nc,x,ix,jx,y,g,ww,ju,vq,cl,ne,nx,nl   2132 
     *am,flmin,  ulam,thr,intr,maxit,xv,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,
     *alm,nlp,jerr)
      goto 14911                                                           2133
14871 continue                                                             2133
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs)                2134
      if(isd .le. 0)goto 14931                                             2134
14940 do 14941 j=1,ni                                                      2134
      cl(:,j)=cl(:,j)*xs(j)                                                2134
14941 continue                                                             2134
14942 continue                                                             2134
14931 continue                                                             2135
      call sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,ww,ju,vq,cl,ne,nx,nlam,f   2138 
     *lmin,  ulam,thr,isd,intr,maxit,kopt,xm,xs,lmu,a0,ca,  ia,nin,dev0,
     *dev,alm,nlp,jerr)
14911 continue                                                             2139
14821 continue                                                             2139
      if(jerr.gt.0) return                                                 2139
      dev0=2.0*sw*dev0                                                     2140
14950 do 14951 k=1,lmu                                                     2140
      nk=nin(k)                                                            2141
14960 do 14961 ic=1,nc                                                     2141
      if(isd .le. 0)goto 14981                                             2141
14990 do 14991 l=1,nk                                                      2141
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                      2141
14991 continue                                                             2141
14992 continue                                                             2141
14981 continue                                                             2142
      if(intr .ne. 0)goto 15011                                            2142
      a0(ic,k)=0.0                                                         2142
      goto 15021                                                           2143
15011 continue                                                             2143
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))            2143
15021 continue                                                             2144
15001 continue                                                             2144
14961 continue                                                             2145
14962 continue                                                             2145
14951 continue                                                             2146
14952 continue                                                             2146
      deallocate(ww,ju,vq,xm,xs)                                           2146
      if(kopt.eq.2) deallocate(xv)                                         2147
      return                                                               2148
      end                                                                  2149
      subroutine multsplstandard2(no,ni,x,ix,jx,w,ju,isd,intr,xm,xs,xv)    2150
      implicit double precision(a-h,o-z)                                   2151
      double precision x(*),w(no),xm(ni),xs(ni),xv(ni)                     2151
      integer ix(*),jx(*),ju(ni)                                           2152
      if(intr .ne. 0)goto 15041                                            2153
15050 do 15051 j=1,ni                                                      2153
      if(ju(j).eq.0)goto 15051                                             2153
      xm(j)=0.0                                                            2153
      jb=ix(j)                                                             2153
      je=ix(j+1)-1                                                         2154
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)                          2155
      if(isd .eq. 0)goto 15071                                             2155
      xbq=dot_product(w(jx(jb:je)),x(jb:je))**2                            2155
      vc=xv(j)-xbq                                                         2156
      xs(j)=sqrt(vc)                                                       2156
      xv(j)=1.0+xbq/vc                                                     2157
      goto 15081                                                           2158
15071 continue                                                             2158
      xs(j)=1.0                                                            2158
15081 continue                                                             2159
15061 continue                                                             2159
15051 continue                                                             2160
15052 continue                                                             2160
      return                                                               2161
15041 continue                                                             2162
15090 do 15091 j=1,ni                                                      2162
      if(ju(j).eq.0)goto 15091                                             2162
      jb=ix(j)                                                             2162
      je=ix(j+1)-1                                                         2163
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2164
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 2165
      if(isd .le. 0)goto 15111                                             2165
      xs(j)=sqrt(xv(j))                                                    2165
      xv(j)=1.0                                                            2165
15111 continue                                                             2166
15091 continue                                                             2167
15092 continue                                                             2167
      if(isd.eq.0) xs=1.0                                                  2168
      return                                                               2169
      end                                                                  2170
      subroutine splstandard2(no,ni,x,ix,jx,w,ju,isd,intr,xm,xs)           2171
      implicit double precision(a-h,o-z)                                   2172
      double precision x(*),w(no),xm(ni),xs(ni)                            2172
      integer ix(*),jx(*),ju(ni)                                           2173
      if(intr .ne. 0)goto 15131                                            2174
15140 do 15141 j=1,ni                                                      2174
      if(ju(j).eq.0)goto 15141                                             2174
      xm(j)=0.0                                                            2174
      jb=ix(j)                                                             2174
      je=ix(j+1)-1                                                         2175
      if(isd .eq. 0)goto 15161                                             2176
      vc=dot_product(w(jx(jb:je)),x(jb:je)**2)  -dot_product(w(jx(jb:je)   2178 
     *),x(jb:je))**2
      xs(j)=sqrt(vc)                                                       2179
      goto 15171                                                           2180
15161 continue                                                             2180
      xs(j)=1.0                                                            2180
15171 continue                                                             2181
15151 continue                                                             2181
15141 continue                                                             2182
15142 continue                                                             2182
      return                                                               2183
15131 continue                                                             2184
15180 do 15181 j=1,ni                                                      2184
      if(ju(j).eq.0)goto 15181                                             2184
      jb=ix(j)                                                             2184
      je=ix(j+1)-1                                                         2185
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2186
      if(isd.ne.0) xs(j)=sqrt(dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j   2187 
     *)**2)
15181 continue                                                             2188
15182 continue                                                             2188
      if(isd.eq.0) xs=1.0                                                  2189
      return                                                               2190
      end                                                                  2191
      subroutine sprlognet2n (parm,no,ni,x,ix,jx,y,g,w,ju,vp,cl,ne,nx,nl   2194 
     *am,  flmin,ulam,shri,isd,intr,maxit,kopt,xb,xs,  lmu,a0,a,m,kin,de
     *v0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   2195
      double precision x(*),y(no),g(no),w(no),vp(ni),ulam(nlam),cl(2,ni)   2196
      double precision a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)             2197
      double precision xb(ni),xs(ni)                                       2197
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           2198
      double precision, dimension (:), allocatable :: xm,b,bs,v,r               
      double precision, dimension (:), allocatable :: sc,xv,q,ga                
      integer, dimension (:), allocatable :: mm,ixx                             
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               2204
      allocate(b(0:ni),stat=jerr)                                          2205
      if(jerr.ne.0) return                                                 2206
      allocate(xm(0:ni),stat=jerr)                                         2207
      if(jerr.ne.0) return                                                 2208
      allocate(xv(1:ni),stat=jerr)                                         2209
      if(jerr.ne.0) return                                                 2210
      allocate(bs(0:ni),stat=jerr)                                         2211
      if(jerr.ne.0) return                                                 2212
      allocate(ga(1:ni),stat=jerr)                                         2213
      if(jerr.ne.0) return                                                 2214
      allocate(mm(1:ni),stat=jerr)                                         2215
      if(jerr.ne.0) return                                                 2216
      allocate(ixx(1:ni),stat=jerr)                                        2217
      if(jerr.ne.0) return                                                 2218
      allocate(q(1:no),stat=jerr)                                          2219
      if(jerr.ne.0) return                                                 2220
      allocate(r(1:no),stat=jerr)                                          2221
      if(jerr.ne.0) return                                                 2222
      allocate(v(1:no),stat=jerr)                                          2223
      if(jerr.ne.0) return                                                 2224
      allocate(sc(1:no),stat=jerr)                                         2225
      if(jerr.ne.0) return                                                 2226
      fmax=log(1.0/pmin-1.0)                                               2226
      fmin=-fmax                                                           2226
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                      2227
      bta=parm                                                             2227
      omb=1.0-bta                                                          2228
      q0=dot_product(w,y)                                                  2228
      if(q0 .gt. pmin)goto 15201                                           2228
      jerr=8001                                                            2228
      return                                                               2228
15201 continue                                                             2229
      if(q0 .lt. 1.0-pmin)goto 15221                                       2229
      jerr=9001                                                            2229
      return                                                               2229
15221 continue                                                             2230
      if(intr.eq.0) q0=0.5                                                 2230
      bz=0.0                                                               2230
      if(intr.ne.0) bz=log(q0/(1.0-q0))                                    2231
      if(nonzero(no,g) .ne. 0)goto 15241                                   2231
      vi=q0*(1.0-q0)                                                       2231
      b(0)=bz                                                              2231
      v=vi*w                                                               2232
      r=w*(y-q0)                                                           2232
      q=q0                                                                 2232
      xm(0)=vi                                                             2232
      dev1=-(bz*q0+log(1.0-q0))                                            2233
      goto 15251                                                           2234
15241 continue                                                             2234
      b(0)=0.0                                                             2235
      if(intr .eq. 0)goto 15271                                            2235
      b(0)=azero(no,y,g,w,jerr)                                            2235
      if(jerr.ne.0) return                                                 2235
15271 continue                                                             2236
      q=1.0/(1.0+exp(-b(0)-g))                                             2236
      v=w*q*(1.0-q)                                                        2236
      r=w*(y-q)                                                            2236
      xm(0)=sum(v)                                                         2237
      dev1=-(b(0)*q0+dot_product(w,y*g+log(1.0-q)))                        2238
15251 continue                                                             2239
15231 continue                                                             2239
      if(kopt .le. 0)goto 15291                                            2240
      if(isd .le. 0 .or. intr .eq. 0)goto 15311                            2240
      xv=0.25                                                              2240
      goto 15321                                                           2241
15311 continue                                                             2242
15330 do 15331 j=1,ni                                                      2242
      if(ju(j).eq.0)goto 15331                                             2242
      jb=ix(j)                                                             2242
      je=ix(j+1)-1                                                         2243
      xv(j)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)          2244
15331 continue                                                             2245
15332 continue                                                             2245
15321 continue                                                             2246
15301 continue                                                             2246
15291 continue                                                             2247
      b(1:ni)=0.0                                                          2247
      dev0=dev1                                                            2248
15340 do 15341 i=1,no                                                      2248
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                        2249
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))              2250
15341 continue                                                             2252
15342 continue                                                             2252
      alf=1.0                                                              2254
      if(flmin .ge. 1.0)goto 15361                                         2254
      eqs=max(eps,flmin)                                                   2254
      alf=eqs**(1.0/(nlam-1))                                              2254
15361 continue                                                             2255
      m=0                                                                  2255
      mm=0                                                                 2255
      nin=0                                                                2255
      o=0.0                                                                2255
      svr=o                                                                2255
      mnl=min(mnlam,nlam)                                                  2255
      bs=0.0                                                               2255
      nlp=0                                                                2255
      nin=nlp                                                              2256
      shr=shri*dev0                                                        2256
      al=0.0                                                               2256
      ixx=0                                                                2257
15370 do 15371 j=1,ni                                                      2257
      if(ju(j).eq.0)goto 15371                                             2258
      jb=ix(j)                                                             2258
      je=ix(j+1)-1                                                         2258
      jn=ix(j+1)-ix(j)                                                     2259
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 2260
      gj=dot_product(sc(1:jn),x(jb:je))                                    2261
      ga(j)=abs((gj-svr*xb(j))/xs(j))                                      2262
15371 continue                                                             2263
15372 continue                                                             2263
15380 do 15381 ilm=1,nlam                                                  2263
      al0=al                                                               2264
      if(flmin .lt. 1.0)goto 15401                                         2264
      al=ulam(ilm)                                                         2264
      goto 15391                                                           2265
15401 if(ilm .le. 2)goto 15411                                             2265
      al=al*alf                                                            2265
      goto 15391                                                           2266
15411 if(ilm .ne. 1)goto 15421                                             2266
      al=big                                                               2266
      goto 15431                                                           2267
15421 continue                                                             2267
      al0=0.0                                                              2268
15440 do 15441 j=1,ni                                                      2268
      if(ju(j).eq.0)goto 15441                                             2268
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2268
15441 continue                                                             2269
15442 continue                                                             2269
      al0=al0/max(bta,1.0d-3)                                              2269
      al=alf*al0                                                           2270
15431 continue                                                             2271
15391 continue                                                             2271
      al2=al*omb                                                           2271
      al1=al*bta                                                           2271
      tlam=bta*(2.0*al-al0)                                                2272
15450 do 15451 k=1,ni                                                      2272
      if(ixx(k).eq.1)goto 15451                                            2272
      if(ju(k).eq.0)goto 15451                                             2273
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     2274
15451 continue                                                             2275
15452 continue                                                             2275
11020 continue                                                             2276
15460 continue                                                             2276
15461 continue                                                             2276
      bs(0)=b(0)                                                           2276
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                2277
15470 do 15471 j=1,ni                                                      2277
      if(ixx(j).eq.0)goto 15471                                            2278
      jb=ix(j)                                                             2278
      je=ix(j+1)-1                                                         2278
      jn=ix(j+1)-ix(j)                                                     2279
      sc(1:jn)=v(jx(jb:je))                                                2280
      xm(j)=dot_product(sc(1:jn),x(jb:je))                                 2281
      if(kopt .ne. 0)goto 15491                                            2282
      xv(j)=dot_product(sc(1:jn),x(jb:je)**2)                              2283
      xv(j)=(xv(j)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2                2284
15491 continue                                                             2285
15471 continue                                                             2286
15472 continue                                                             2286
15500 continue                                                             2286
15501 continue                                                             2286
      nlp=nlp+1                                                            2286
      dlx=0.0                                                              2287
15510 do 15511 k=1,ni                                                      2287
      if(ixx(k).eq.0)goto 15511                                            2288
      jb=ix(k)                                                             2288
      je=ix(k+1)-1                                                         2288
      jn=ix(k+1)-ix(k)                                                     2288
      bk=b(k)                                                              2289
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 2290
      gk=dot_product(sc(1:jn),x(jb:je))                                    2291
      gk=(gk-svr*xb(k))/xs(k)                                              2292
      u=gk+xv(k)*b(k)                                                      2292
      au=abs(u)-vp(k)*al1                                                  2293
      if(au .gt. 0.0)goto 15531                                            2293
      b(k)=0.0                                                             2293
      goto 15541                                                           2294
15531 continue                                                             2295
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))          2296
15541 continue                                                             2297
15521 continue                                                             2297
      d=b(k)-bk                                                            2297
      if(abs(d).le.0.0)goto 15511                                          2297
      dlx=max(dlx,xv(k)*d**2)                                              2298
      if(mm(k) .ne. 0)goto 15561                                           2298
      nin=nin+1                                                            2298
      if(nin.gt.nx)goto 15512                                              2299
      mm(k)=nin                                                            2299
      m(nin)=k                                                             2299
      sc(1:jn)=v(jx(jb:je))                                                2300
      xm(k)=dot_product(sc(1:jn),x(jb:je))                                 2301
15561 continue                                                             2302
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              2303
      o=o+d*(xb(k)/xs(k))                                                  2304
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  2305
15511 continue                                                             2306
15512 continue                                                             2306
      if(nin.gt.nx)goto 15502                                              2307
      d=0.0                                                                2307
      if(intr.ne.0) d=svr/xm(0)                                            2308
      if(d .eq. 0.0)goto 15581                                             2308
      b(0)=b(0)+d                                                          2308
      dlx=max(dlx,xm(0)*d**2)                                              2308
      r=r-d*v                                                              2309
      svr=svr-d*xm(0)                                                      2310
15581 continue                                                             2311
      if(dlx.lt.shr)goto 15502                                             2312
      if(nlp .le. maxit)goto 15601                                         2312
      jerr=-ilm                                                            2312
      return                                                               2312
15601 continue                                                             2313
15610 continue                                                             2313
15611 continue                                                             2313
      nlp=nlp+1                                                            2313
      dlx=0.0                                                              2314
15620 do 15621 l=1,nin                                                     2314
      k=m(l)                                                               2314
      jb=ix(k)                                                             2314
      je=ix(k+1)-1                                                         2315
      jn=ix(k+1)-ix(k)                                                     2315
      bk=b(k)                                                              2316
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 2317
      gk=dot_product(sc(1:jn),x(jb:je))                                    2318
      gk=(gk-svr*xb(k))/xs(k)                                              2319
      u=gk+xv(k)*b(k)                                                      2319
      au=abs(u)-vp(k)*al1                                                  2320
      if(au .gt. 0.0)goto 15641                                            2320
      b(k)=0.0                                                             2320
      goto 15651                                                           2321
15641 continue                                                             2322
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))          2323
15651 continue                                                             2324
15631 continue                                                             2324
      d=b(k)-bk                                                            2324
      if(abs(d).le.0.0)goto 15621                                          2324
      dlx=max(dlx,xv(k)*d**2)                                              2325
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              2326
      o=o+d*(xb(k)/xs(k))                                                  2327
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  2328
15621 continue                                                             2329
15622 continue                                                             2329
      d=0.0                                                                2329
      if(intr.ne.0) d=svr/xm(0)                                            2330
      if(d .eq. 0.0)goto 15671                                             2330
      b(0)=b(0)+d                                                          2330
      dlx=max(dlx,xm(0)*d**2)                                              2330
      r=r-d*v                                                              2331
      svr=svr-d*xm(0)                                                      2332
15671 continue                                                             2333
      if(dlx.lt.shr)goto 15612                                             2334
      if(nlp .le. maxit)goto 15691                                         2334
      jerr=-ilm                                                            2334
      return                                                               2334
15691 continue                                                             2335
      goto 15611                                                           2336
15612 continue                                                             2336
      goto 15501                                                           2337
15502 continue                                                             2337
      if(nin.gt.nx)goto 15462                                              2338
      sc=b(0)                                                              2338
      b0=0.0                                                               2339
15700 do 15701 j=1,nin                                                     2339
      l=m(j)                                                               2339
      jb=ix(l)                                                             2339
      je=ix(l+1)-1                                                         2340
      sc(jx(jb:je))=sc(jx(jb:je))+b(l)*x(jb:je)/xs(l)                      2341
      b0=b0-b(l)*xb(l)/xs(l)                                               2342
15701 continue                                                             2343
15702 continue                                                             2343
      sc=sc+b0                                                             2344
15710 do 15711 i=1,no                                                      2344
      fi=sc(i)+g(i)                                                        2345
      if(fi .ge. fmin)goto 15731                                           2345
      q(i)=0.0                                                             2345
      goto 15721                                                           2345
15731 if(fi .le. fmax)goto 15741                                           2345
      q(i)=1.0                                                             2345
      goto 15751                                                           2346
15741 continue                                                             2346
      q(i)=1.0/(1.0+exp(-fi))                                              2346
15751 continue                                                             2347
15721 continue                                                             2347
15711 continue                                                             2348
15712 continue                                                             2348
      v=w*q*(1.0-q)                                                        2348
      xm(0)=sum(v)                                                         2348
      if(xm(0).lt.vmin)goto 15462                                          2349
      r=w*(y-q)                                                            2349
      svr=sum(r)                                                           2349
      o=0.0                                                                2350
      if(xm(0)*(b(0)-bs(0))**2 .ge. shr)goto 15771                         2350
      kx=0                                                                 2351
15780 do 15781 j=1,nin                                                     2351
      k=m(j)                                                               2352
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 15781                           2352
      kx=1                                                                 2352
      goto 15782                                                           2353
15781 continue                                                             2354
15782 continue                                                             2354
      if(kx .ne. 0)goto 15801                                              2355
15810 do 15811 j=1,ni                                                      2355
      if(ixx(j).eq.1)goto 15811                                            2355
      if(ju(j).eq.0)goto 15811                                             2356
      jb=ix(j)                                                             2356
      je=ix(j+1)-1                                                         2356
      jn=ix(j+1)-ix(j)                                                     2357
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 2358
      gj=dot_product(sc(1:jn),x(jb:je))                                    2359
      ga(j)=abs((gj-svr*xb(j))/xs(j))                                      2360
      if(ga(j) .le. al1*vp(j))goto 15831                                   2360
      ixx(j)=1                                                             2360
      kx=1                                                                 2360
15831 continue                                                             2361
15811 continue                                                             2362
15812 continue                                                             2362
      if(kx.eq.1) go to 11020                                              2363
      goto 15462                                                           2364
15801 continue                                                             2365
15771 continue                                                             2366
      goto 15461                                                           2367
15462 continue                                                             2367
      if(nin .le. nx)goto 15851                                            2367
      jerr=-10000-ilm                                                      2367
      goto 15382                                                           2367
15851 continue                                                             2368
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                2368
      kin(ilm)=nin                                                         2369
      a0(ilm)=b(0)                                                         2369
      alm(ilm)=al                                                          2369
      lmu=ilm                                                              2370
      devi=dev2(no,w,y,q,pmin)                                             2371
      dev(ilm)=(dev1-devi)/dev0                                            2372
      if(ilm.lt.mnl)goto 15381                                             2372
      if(flmin.ge.1.0)goto 15381                                           2373
      me=0                                                                 2373
15860 do 15861 j=1,nin                                                     2373
      if(a(j,ilm).ne.0.0) me=me+1                                          2373
15861 continue                                                             2373
15862 continue                                                             2373
      if(me.gt.ne)goto 15382                                               2374
      if(dev(ilm).gt.devmax)goto 15382                                     2374
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 15382                             2375
      if(xm(0).lt.vmin)goto 15382                                          2376
15381 continue                                                             2377
15382 continue                                                             2377
      g=log(q/(1.0-q))                                                     2378
      deallocate(xm,b,bs,v,r,sc,xv,q,mm,ga,ixx)                            2379
      return                                                               2380
      end                                                                  2381
      subroutine sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,w,ju,vp,cl,ne,nx,n   2383 
     *lam,flmin,  ulam,shri,isd,intr,maxit,kopt,xb,xs,lmu,a0,a,m,kin,dev
     *0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   2384
      double precision x(*),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam),xb   2385 
     *(ni),xs(ni)
      double precision a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(   2386 
     *2,ni)
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           2387
      double precision, dimension (:,:), allocatable :: q                       
      double precision, dimension (:), allocatable :: sxp,sxpl                  
      double precision, dimension (:), allocatable :: sc,xm,v,r,ga              
      double precision, dimension (:,:), allocatable :: b,bs,xv                 
      integer, dimension (:), allocatable :: mm,is,iy                           
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      if(jerr.ne.0) return                                                      
      allocate(xv(1:ni,1:nc),stat=jerr)                                         
      if(jerr.ne.0) return 					                                                
      allocate(bs(0:ni,1:nc),stat=jerr)                                         
      if(jerr.ne.0) return					                                                 
      allocate(q(1:no,1:nc),stat=jerr)                                          
      if(jerr.ne.0) return					                                                 
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               2402
      exmn=-exmx                                                           2403
      allocate(xm(0:ni),stat=jerr)                                         2404
      if(jerr.ne.0) return                                                 2405
      allocate(r(1:no),stat=jerr)                                          2406
      if(jerr.ne.0) return                                                 2407
      allocate(v(1:no),stat=jerr)                                          2408
      if(jerr.ne.0) return                                                 2409
      allocate(mm(1:ni),stat=jerr)                                         2410
      if(jerr.ne.0) return                                                 2411
      allocate(ga(1:ni),stat=jerr)                                         2412
      if(jerr.ne.0) return                                                 2413
      allocate(iy(1:ni),stat=jerr)                                         2414
      if(jerr.ne.0) return                                                 2415
      allocate(is(1:max(nc,ni)),stat=jerr)                                 2416
      if(jerr.ne.0) return                                                 2417
      allocate(sxp(1:no),stat=jerr)                                        2418
      if(jerr.ne.0) return                                                 2419
      allocate(sxpl(1:no),stat=jerr)                                       2420
      if(jerr.ne.0) return                                                 2421
      allocate(sc(1:no),stat=jerr)                                         2422
      if(jerr.ne.0) return                                                 2423
      pmax=1.0-pmin                                                        2423
      emin=pmin/pmax                                                       2423
      emax=1.0/emin                                                        2424
      pfm=(1.0+pmin)*pmin                                                  2424
      pfx=(1.0-pmin)*pmax                                                  2424
      vmin=pfm*pmax                                                        2425
      bta=parm                                                             2425
      omb=1.0-bta                                                          2425
      dev1=0.0                                                             2425
      dev0=0.0                                                             2426
15870 do 15871 ic=1,nc                                                     2426
      q0=dot_product(w,y(:,ic))                                            2427
      if(q0 .gt. pmin)goto 15891                                           2427
      jerr =8000+ic                                                        2427
      return                                                               2427
15891 continue                                                             2428
      if(q0 .lt. 1.0-pmin)goto 15911                                       2428
      jerr =9000+ic                                                        2428
      return                                                               2428
15911 continue                                                             2429
      if(intr.eq.0) q0=1.0/nc                                              2430
      b(1:ni,ic)=0.0                                                       2430
      b(0,ic)=0.0                                                          2431
      if(intr .eq. 0)goto 15931                                            2431
      b(0,ic)=log(q0)                                                      2431
      dev1=dev1-q0*b(0,ic)                                                 2431
15931 continue                                                             2432
15871 continue                                                             2433
15872 continue                                                             2433
      if(intr.eq.0) dev1=log(float(nc))                                    2433
      iy=0                                                                 2433
      al=0.0                                                               2434
      if(nonzero(no*nc,g) .ne. 0)goto 15951                                2435
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         2435
      sxp=0.0                                                              2436
15960 do 15961 ic=1,nc                                                     2436
      q(:,ic)=exp(b(0,ic))                                                 2436
      sxp=sxp+q(:,ic)                                                      2436
15961 continue                                                             2437
15962 continue                                                             2437
      goto 15971                                                           2438
15951 continue                                                             2438
15980 do 15981 i=1,no                                                      2438
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         2438
15981 continue                                                             2438
15982 continue                                                             2438
      sxp=0.0                                                              2439
      if(intr .ne. 0)goto 16001                                            2439
      b(0,:)=0.0                                                           2439
      goto 16011                                                           2440
16001 continue                                                             2440
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 2440
      if(jerr.ne.0) return                                                 2440
16011 continue                                                             2441
15991 continue                                                             2441
      dev1=0.0                                                             2442
16020 do 16021 ic=1,nc                                                     2442
      q(:,ic)=b(0,ic)+g(:,ic)                                              2443
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             2444
      q(:,ic)=exp(q(:,ic))                                                 2444
      sxp=sxp+q(:,ic)                                                      2445
16021 continue                                                             2446
16022 continue                                                             2446
      sxpl=w*log(sxp)                                                      2446
16030 do 16031 ic=1,nc                                                     2446
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  2446
16031 continue                                                             2447
16032 continue                                                             2447
15971 continue                                                             2448
15941 continue                                                             2448
16040 do 16041 ic=1,nc                                                     2448
16050 do 16051 i=1,no                                                      2448
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               2448
16051 continue                                                             2448
16052 continue                                                             2448
16041 continue                                                             2449
16042 continue                                                             2449
      dev0=dev0+dev1                                                       2450
      if(kopt .le. 0)goto 16071                                            2451
      if(isd .le. 0 .or. intr .eq. 0)goto 16091                            2451
      xv=0.25                                                              2451
      goto 16101                                                           2452
16091 continue                                                             2453
16110 do 16111 j=1,ni                                                      2453
      if(ju(j).eq.0)goto 16111                                             2453
      jb=ix(j)                                                             2453
      je=ix(j+1)-1                                                         2454
      xv(j,:)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)        2455
16111 continue                                                             2456
16112 continue                                                             2456
16101 continue                                                             2457
16081 continue                                                             2457
16071 continue                                                             2459
      alf=1.0                                                              2461
      if(flmin .ge. 1.0)goto 16131                                         2461
      eqs=max(eps,flmin)                                                   2461
      alf=eqs**(1.0/(nlam-1))                                              2461
16131 continue                                                             2462
      m=0                                                                  2462
      mm=0                                                                 2462
      nin=0                                                                2462
      nlp=0                                                                2462
      mnl=min(mnlam,nlam)                                                  2462
      bs=0.0                                                               2462
      svr=0.0                                                              2462
      o=0.0                                                                2463
      shr=shri*dev0                                                        2463
      ga=0.0                                                               2464
16140 do 16141 ic=1,nc                                                     2464
      v=q(:,ic)/sxp                                                        2464
      r=w*(y(:,ic)-v)                                                      2464
      v=w*v*(1.0-v)                                                        2465
16150 do 16151 j=1,ni                                                      2465
      if(ju(j).eq.0)goto 16151                                             2466
      jb=ix(j)                                                             2466
      je=ix(j+1)-1                                                         2466
      jn=ix(j+1)-ix(j)                                                     2467
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2468
      gj=dot_product(sc(1:jn),x(jb:je))                                    2469
      ga(j)=max(ga(j),abs(gj-svr*xb(j))/xs(j))                             2470
16151 continue                                                             2471
16152 continue                                                             2471
16141 continue                                                             2472
16142 continue                                                             2472
16160 do 16161 ilm=1,nlam                                                  2472
      al0=al                                                               2473
      if(flmin .lt. 1.0)goto 16181                                         2473
      al=ulam(ilm)                                                         2473
      goto 16171                                                           2474
16181 if(ilm .le. 2)goto 16191                                             2474
      al=al*alf                                                            2474
      goto 16171                                                           2475
16191 if(ilm .ne. 1)goto 16201                                             2475
      al=big                                                               2475
      goto 16211                                                           2476
16201 continue                                                             2476
      al0=0.0                                                              2477
16220 do 16221 j=1,ni                                                      2477
      if(ju(j).eq.0)goto 16221                                             2477
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2477
16221 continue                                                             2478
16222 continue                                                             2478
      al0=al0/max(bta,1.0d-3)                                              2478
      al=alf*al0                                                           2479
16211 continue                                                             2480
16171 continue                                                             2480
      al2=al*omb                                                           2480
      al1=al*bta                                                           2480
      tlam=bta*(2.0*al-al0)                                                2481
16230 do 16231 k=1,ni                                                      2481
      if(iy(k).eq.1)goto 16231                                             2481
      if(ju(k).eq.0)goto 16231                                             2482
      if(ga(k).gt.tlam*vp(k)) iy(k)=1                                      2483
16231 continue                                                             2484
16232 continue                                                             2484
11020 continue                                                             2485
16240 continue                                                             2485
16241 continue                                                             2485
      ixx=0                                                                2485
      jxx=ixx                                                              2485
      ig=0                                                                 2486
16250 do 16251 ic=1,nc                                                     2486
      bs(0,ic)=b(0,ic)                                                     2487
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          2488
      xm(0)=0.0                                                            2488
      svr=0.0                                                              2488
      o=0.0                                                                2489
16260 do 16261 i=1,no                                                      2489
      pic=q(i,ic)/sxp(i)                                                   2490
      if(pic .ge. pfm)goto 16281                                           2490
      pic=0.0                                                              2490
      v(i)=0.0                                                             2490
      goto 16271                                                           2491
16281 if(pic .le. pfx)goto 16291                                           2491
      pic=1.0                                                              2491
      v(i)=0.0                                                             2491
      goto 16301                                                           2492
16291 continue                                                             2492
      v(i)=w(i)*pic*(1.0-pic)                                              2492
      xm(0)=xm(0)+v(i)                                                     2492
16301 continue                                                             2493
16271 continue                                                             2493
      r(i)=w(i)*(y(i,ic)-pic)                                              2493
      svr=svr+r(i)                                                         2494
16261 continue                                                             2495
16262 continue                                                             2495
      if(xm(0).le.vmin)goto 16251                                          2495
      ig=1                                                                 2496
16310 do 16311 j=1,ni                                                      2496
      if(iy(j).eq.0)goto 16311                                             2497
      jb=ix(j)                                                             2497
      je=ix(j+1)-1                                                         2498
      xm(j)=dot_product(v(jx(jb:je)),x(jb:je))                             2499
      if(kopt .ne. 0)goto 16331                                            2500
      xv(j,ic)=dot_product(v(jx(jb:je)),x(jb:je)**2)                       2501
      xv(j,ic)=(xv(j,ic)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2          2502
16331 continue                                                             2503
16311 continue                                                             2504
16312 continue                                                             2504
16340 continue                                                             2504
16341 continue                                                             2504
      nlp=nlp+1                                                            2504
      dlx=0.0                                                              2505
16350 do 16351 k=1,ni                                                      2505
      if(iy(k).eq.0)goto 16351                                             2506
      jb=ix(k)                                                             2506
      je=ix(k+1)-1                                                         2506
      jn=ix(k+1)-ix(k)                                                     2506
      bk=b(k,ic)                                                           2507
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2508
      gk=dot_product(sc(1:jn),x(jb:je))                                    2509
      gk=(gk-svr*xb(k))/xs(k)                                              2510
      u=gk+xv(k,ic)*b(k,ic)                                                2510
      au=abs(u)-vp(k)*al1                                                  2511
      if(au .gt. 0.0)goto 16371                                            2511
      b(k,ic)=0.0                                                          2511
      goto 16381                                                           2512
16371 continue                                                             2513
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))   2515 
     *)
16381 continue                                                             2516
16361 continue                                                             2516
      d=b(k,ic)-bk                                                         2516
      if(abs(d).le.0.0)goto 16351                                          2517
      dlx=max(dlx,xv(k,ic)*d**2)                                           2518
      if(mm(k) .ne. 0)goto 16401                                           2518
      nin=nin+1                                                            2519
      if(nin .le. nx)goto 16421                                            2519
      jxx=1                                                                2519
      goto 16352                                                           2519
16421 continue                                                             2520
      mm(k)=nin                                                            2520
      m(nin)=k                                                             2521
      xm(k)=dot_product(v(jx(jb:je)),x(jb:je))                             2522
16401 continue                                                             2523
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              2524
      o=o+d*(xb(k)/xs(k))                                                  2525
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  2526
16351 continue                                                             2527
16352 continue                                                             2527
      if(jxx.gt.0)goto 16342                                               2528
      d=0.0                                                                2528
      if(intr.ne.0) d=svr/xm(0)                                            2529
      if(d .eq. 0.0)goto 16441                                             2529
      b(0,ic)=b(0,ic)+d                                                    2529
      dlx=max(dlx,xm(0)*d**2)                                              2530
      r=r-d*v                                                              2530
      svr=svr-d*xm(0)                                                      2531
16441 continue                                                             2532
      if(dlx.lt.shr)goto 16342                                             2532
      if(nlp .le. maxit)goto 16461                                         2532
      jerr=-ilm                                                            2532
      return                                                               2532
16461 continue                                                             2533
16470 continue                                                             2533
16471 continue                                                             2533
      nlp=nlp+1                                                            2533
      dlx=0.0                                                              2534
16480 do 16481 l=1,nin                                                     2534
      k=m(l)                                                               2534
      jb=ix(k)                                                             2534
      je=ix(k+1)-1                                                         2535
      jn=ix(k+1)-ix(k)                                                     2535
      bk=b(k,ic)                                                           2536
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2537
      gk=dot_product(sc(1:jn),x(jb:je))                                    2538
      gk=(gk-svr*xb(k))/xs(k)                                              2539
      u=gk+xv(k,ic)*b(k,ic)                                                2539
      au=abs(u)-vp(k)*al1                                                  2540
      if(au .gt. 0.0)goto 16501                                            2540
      b(k,ic)=0.0                                                          2540
      goto 16511                                                           2541
16501 continue                                                             2542
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))   2544 
     *)
16511 continue                                                             2545
16491 continue                                                             2545
      d=b(k,ic)-bk                                                         2545
      if(abs(d).le.0.0)goto 16481                                          2546
      dlx=max(dlx,xv(k,ic)*d**2)                                           2547
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              2548
      o=o+d*(xb(k)/xs(k))                                                  2549
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  2550
16481 continue                                                             2551
16482 continue                                                             2551
      d=0.0                                                                2551
      if(intr.ne.0) d=svr/xm(0)                                            2552
      if(d .eq. 0.0)goto 16531                                             2552
      b(0,ic)=b(0,ic)+d                                                    2552
      dlx=max(dlx,xm(0)*d**2)                                              2553
      r=r-d*v                                                              2553
      svr=svr-d*xm(0)                                                      2554
16531 continue                                                             2555
      if(dlx.lt.shr)goto 16472                                             2555
      if(nlp .le. maxit)goto 16551                                         2555
      jerr=-ilm                                                            2555
      return                                                               2555
16551 continue                                                             2556
      goto 16471                                                           2557
16472 continue                                                             2557
      goto 16341                                                           2558
16342 continue                                                             2558
      if(jxx.gt.0)goto 16252                                               2559
      if(xm(0)*(b(0,ic)-bs(0,ic))**2.gt.shr) ixx=1                         2560
      if(ixx .ne. 0)goto 16571                                             2561
16580 do 16581 j=1,nin                                                     2561
      k=m(j)                                                               2562
      if(xv(k,ic)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 16601                2562
      ixx=1                                                                2562
      goto 16582                                                           2562
16601 continue                                                             2563
16581 continue                                                             2564
16582 continue                                                             2564
16571 continue                                                             2565
      sc=b(0,ic)+g(:,ic)                                                   2565
      b0=0.0                                                               2566
16610 do 16611 j=1,nin                                                     2566
      l=m(j)                                                               2566
      jb=ix(l)                                                             2566
      je=ix(l+1)-1                                                         2567
      sc(jx(jb:je))=sc(jx(jb:je))+b(l,ic)*x(jb:je)/xs(l)                   2568
      b0=b0-b(l,ic)*xb(l)/xs(l)                                            2569
16611 continue                                                             2570
16612 continue                                                             2570
      sc=min(max(exmn,sc+b0),exmx)                                         2571
      sxp=sxp-q(:,ic)                                                      2572
      q(:,ic)=min(max(emin*sxp,exp(sc)),emax*sxp)                          2573
      sxp=sxp+q(:,ic)                                                      2574
16251 continue                                                             2575
16252 continue                                                             2575
      s=-sum(b(0,:))/nc                                                    2575
      b(0,:)=b(0,:)+s                                                      2575
      sc=s                                                                 2575
      b0=0.0                                                               2576
16620 do 16621 j=1,nin                                                     2576
      l=m(j)                                                               2577
      if(vp(l) .gt. 0.0)goto 16641                                         2577
      s=sum(b(l,:))/nc                                                     2577
      goto 16651                                                           2578
16641 continue                                                             2578
      s=elc(parm,nc,cl(:,l),b(l,:),is)                                     2578
16651 continue                                                             2579
16631 continue                                                             2579
      b(l,:)=b(l,:)-s                                                      2580
      jb=ix(l)                                                             2580
      je=ix(l+1)-1                                                         2581
      sc(jx(jb:je))=sc(jx(jb:je))-s*x(jb:je)/xs(l)                         2582
      b0=b0+s*xb(l)/xs(l)                                                  2583
16621 continue                                                             2584
16622 continue                                                             2584
      sc=sc+b0                                                             2584
      sc=exp(sc)                                                           2584
      sxp=sxp*sc                                                           2584
16660 do 16661 ic=1,nc                                                     2584
      q(:,ic)=q(:,ic)*sc                                                   2584
16661 continue                                                             2585
16662 continue                                                             2585
      if(jxx.gt.0)goto 16242                                               2585
      if(ig.eq.0)goto 16242                                                2586
      if(ixx .ne. 0)goto 16681                                             2587
16690 do 16691 j=1,ni                                                      2587
      if(iy(j).eq.1)goto 16691                                             2587
      if(ju(j).eq.0)goto 16691                                             2587
      ga(j)=0.0                                                            2587
16691 continue                                                             2588
16692 continue                                                             2588
16700 do 16701 ic=1,nc                                                     2588
      v=q(:,ic)/sxp                                                        2588
      r=w*(y(:,ic)-v)                                                      2588
      v=w*v*(1.0-v)                                                        2589
16710 do 16711 j=1,ni                                                      2589
      if(iy(j).eq.1)goto 16711                                             2589
      if(ju(j).eq.0)goto 16711                                             2590
      jb=ix(j)                                                             2590
      je=ix(j+1)-1                                                         2590
      jn=ix(j+1)-ix(j)                                                     2591
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2592
      gj=dot_product(sc(1:jn),x(jb:je))                                    2593
      ga(j)=max(ga(j),abs(gj-svr*xb(j))/xs(j))                             2594
16711 continue                                                             2595
16712 continue                                                             2595
16701 continue                                                             2596
16702 continue                                                             2596
16720 do 16721 k=1,ni                                                      2596
      if(iy(k).eq.1)goto 16721                                             2596
      if(ju(k).eq.0)goto 16721                                             2597
      if(ga(k) .le. al1*vp(k))goto 16741                                   2597
      iy(k)=1                                                              2597
      ixx=1                                                                2597
16741 continue                                                             2598
16721 continue                                                             2599
16722 continue                                                             2599
      if(ixx.eq.1) go to 11020                                             2600
      goto 16242                                                           2601
16681 continue                                                             2602
      goto 16241                                                           2603
16242 continue                                                             2603
      if(jxx .le. 0)goto 16761                                             2603
      jerr=-10000-ilm                                                      2603
      goto 16162                                                           2603
16761 continue                                                             2603
      devi=0.0                                                             2604
16770 do 16771 ic=1,nc                                                     2605
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          2605
      a0(ic,ilm)=b(0,ic)                                                   2606
16780 do 16781 i=1,no                                                      2606
      if(y(i,ic).le.0.0)goto 16781                                         2607
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           2608
16781 continue                                                             2609
16782 continue                                                             2609
16771 continue                                                             2610
16772 continue                                                             2610
      kin(ilm)=nin                                                         2610
      alm(ilm)=al                                                          2610
      lmu=ilm                                                              2611
      dev(ilm)=(dev1-devi)/dev0                                            2611
      if(ig.eq.0)goto 16162                                                2612
      if(ilm.lt.mnl)goto 16161                                             2612
      if(flmin.ge.1.0)goto 16161                                           2613
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 16162             2614
      if(dev(ilm).gt.devmax)goto 16162                                     2614
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 16162                             2615
16161 continue                                                             2616
16162 continue                                                             2616
      g=log(q)                                                             2616
16790 do 16791 i=1,no                                                      2616
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         2616
16791 continue                                                             2617
16792 continue                                                             2617
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,xm,sc,ga,iy)                      2618
      return                                                               2619
      end                                                                  2620
      subroutine lcmodval(nc,nx,a0,ca,ia,nin,x,ix,jx,n,f)                  2621
      implicit double precision(a-h,o-z)                                   2622
      double precision a0(nc),ca(nx,nc),x(*),f(nc,n)                       2622
      integer ia(*),ix(*),jx(*)                                            2623
16800 do 16801 ic=1,nc                                                     2623
      f(ic,:)=a0(ic)                                                       2623
16801 continue                                                             2624
16802 continue                                                             2624
16810 do 16811 j=1,nin                                                     2624
      k=ia(j)                                                              2624
      kb=ix(k)                                                             2624
      ke=ix(k+1)-1                                                         2625
16820 do 16821 ic=1,nc                                                     2625
      f(ic,jx(kb:ke))=f(ic,jx(kb:ke))+ca(j,ic)*x(kb:ke)                    2625
16821 continue                                                             2626
16822 continue                                                             2626
16811 continue                                                             2627
16812 continue                                                             2627
      return                                                               2628
      end                                                                  2629
      subroutine coxnet (parm,no,ni,x,y,d,g,w,jd,vp,cl,ne,nx,nlam,flmin,   2631 
     *ulam,thr,  maxit,isd,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   2632
      double precision x(no,ni),y(no),d(no),g(no),w(no),vp(ni),ulam(nlam   2633 
     *)
      double precision ca(nx,nlam),dev(nlam),alm(nlam),cl(2,ni)            2634
      integer jd(*),ia(nx),nin(nlam)                                       2635
      double precision, dimension (:), allocatable :: xs,ww,vq                  
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 16841                                    2639
      jerr=10000                                                           2639
      return                                                               2639
16841 continue                                                             2640
      allocate(ww(1:no),stat=jerr)                                         2641
      if(jerr.ne.0) return                                                 2642
      allocate(ju(1:ni),stat=jerr)                                         2643
      if(jerr.ne.0) return                                                 2644
      allocate(vq(1:ni),stat=jerr)                                         2645
      if(jerr.ne.0) return                                                 2646
      if(isd .le. 0)goto 16861                                             2646
      allocate(xs(1:ni),stat=jerr)                                         2646
      if(jerr.ne.0) return                                                 2646
16861 continue                                                             2648
      call chkvars(no,ni,x,ju)                                             2649
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2650
      if(maxval(ju) .gt. 0)goto 16881                                      2650
      jerr=7777                                                            2650
      return                                                               2650
16881 continue                                                             2651
      vq=max(0d0,vp)                                                       2651
      vq=vq*ni/sum(vq)                                                     2652
      ww=max(0d0,w)                                                        2652
      sw=sum(ww)                                                           2653
      if(sw .gt. 0.0)goto 16901                                            2653
      jerr=9999                                                            2653
      return                                                               2653
16901 continue                                                             2653
      ww=ww/sw                                                             2654
      call cstandard(no,ni,x,ww,ju,isd,xs)                                 2655
      if(isd .le. 0)goto 16921                                             2655
16930 do 16931 j=1,ni                                                      2655
      cl(:,j)=cl(:,j)*xs(j)                                                2655
16931 continue                                                             2655
16932 continue                                                             2655
16921 continue                                                             2656
      call coxnet1(parm,no,ni,x,y,d,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam,   2658 
     *thr,  isd,maxit,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 2658
      dev0=2.0*sw*dev0                                                     2659
      if(isd .le. 0)goto 16951                                             2659
16960 do 16961 k=1,lmu                                                     2659
      nk=nin(k)                                                            2659
      ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                                   2659
16961 continue                                                             2659
16962 continue                                                             2659
16951 continue                                                             2660
      deallocate(ww,ju,vq)                                                 2660
      if(isd.gt.0) deallocate(xs)                                          2661
      return                                                               2662
      end                                                                  2663
      subroutine cstandard (no,ni,x,w,ju,isd,xs)                           2664
      implicit double precision(a-h,o-z)                                   2665
      double precision x(no,ni),w(no),xs(ni)                               2665
      integer ju(ni)                                                       2666
16970 do 16971 j=1,ni                                                      2666
      if(ju(j).eq.0)goto 16971                                             2667
      xm=dot_product(w,x(:,j))                                             2667
      x(:,j)=x(:,j)-xm                                                     2668
      if(isd .le. 0)goto 16991                                             2668
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                                 2668
      x(:,j)=x(:,j)/xs(j)                                                  2668
16991 continue                                                             2669
16971 continue                                                             2670
16972 continue                                                             2670
      return                                                               2671
      end                                                                  2672
      subroutine coxnet1(parm,no,ni,x,y,d,g,q,ju,vp,cl,ne,nx,nlam,flmin,   2674 
     *ulam,cthri,  isd,maxit,lmu,ao,m,kin,dev0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   2675
      double precision x(no,ni),y(no),q(no),d(no),g(no),vp(ni),ulam(nlam   2676 
     *)
      double precision ao(nx,nlam),dev(nlam),alm(nlam),cl(2,ni)            2677
      integer ju(ni),m(nx),kin(nlam)                                       2678
      double precision, dimension (:), allocatable :: w,dk,v,xs,wr              
      double precision, dimension (:), allocatable :: a,as,f,dq                 
      double precision, dimension (:), allocatable :: e,uu,ga                   
      integer, dimension (:), allocatable :: jp,kp,mm,ixx                       
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               2685
      sml=sml*100.0                                                        2685
      devmax=devmax*0.99/0.999                                             2686
      allocate(e(1:no),stat=jerr)                                          2687
      if(jerr.ne.0)go to 12320                                             2688
      allocate(uu(1:no),stat=jerr)                                         2689
      if(jerr.ne.0)go to 12320                                             2690
      allocate(f(1:no),stat=jerr)                                          2691
      if(jerr.ne.0)go to 12320                                             2692
      allocate(w(1:no),stat=jerr)                                          2693
      if(jerr.ne.0)go to 12320                                             2694
      allocate(v(1:ni),stat=jerr)                                          2695
      if(jerr.ne.0)go to 12320                                             2696
      allocate(a(1:ni),stat=jerr)                                          2697
      if(jerr.ne.0)go to 12320                                             2698
      allocate(as(1:ni),stat=jerr)                                         2699
      if(jerr.ne.0)go to 12320                                             2700
      allocate(xs(1:ni),stat=jerr)                                         2701
      if(jerr.ne.0)go to 12320                                             2702
      allocate(ga(1:ni),stat=jerr)                                         2703
      if(jerr.ne.0)go to 12320                                             2704
      allocate(ixx(1:ni),stat=jerr)                                        2705
      if(jerr.ne.0)go to 12320                                             2706
      allocate(jp(1:no),stat=jerr)                                         2707
      if(jerr.ne.0)go to 12320                                             2708
      allocate(kp(1:no),stat=jerr)                                         2709
      if(jerr.ne.0)go to 12320                                             2710
      allocate(dk(1:no),stat=jerr)                                         2711
      if(jerr.ne.0)go to 12320                                             2712
      allocate(wr(1:no),stat=jerr)                                         2713
      if(jerr.ne.0)go to 12320                                             2714
      allocate(dq(1:no),stat=jerr)                                         2715
      if(jerr.ne.0)go to 12320                                             2716
      allocate(mm(1:ni),stat=jerr)                                         2717
      if(jerr.ne.0)go to 12320                                             2718
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                               2719
      if(jerr.ne.0) go to 12320                                            2719
      alpha=parm                                                           2720
      oma=1.0-alpha                                                        2720
      nlm=0                                                                2720
      ixx=0                                                                2720
      al=0.0                                                               2721
      dq=d*q                                                               2721
      call died(no,nk,dq,kp,jp,dk)                                         2722
      a=0.0                                                                2722
      f(1)=0.0                                                             2722
      fmax=log(huge(f(1))*0.1)                                             2723
      if(nonzero(no,g) .eq. 0)goto 17011                                   2723
      f=g-dot_product(q,g)                                                 2724
      e=q*exp(sign(min(abs(f),fmax),f))                                    2725
      goto 17021                                                           2726
17011 continue                                                             2726
      f=0.0                                                                2726
      e=q                                                                  2726
17021 continue                                                             2727
17001 continue                                                             2727
      r0=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                                 2728
      rr=-(dot_product(dk(1:nk),log(dk(1:nk)))+r0)                         2728
      dev0=rr                                                              2729
17030 do 17031 i=1,no                                                      2729
      if((y(i) .ge. t0) .and. (q(i) .gt. 0.0))goto 17051                   2729
      w(i)=0.0                                                             2729
      wr(i)=w(i)                                                           2729
17051 continue                                                             2729
17031 continue                                                             2730
17032 continue                                                             2730
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                         2731
      if(jerr.ne.0) go to 12320                                            2733
      alf=1.0                                                              2735
      if(flmin .ge. 1.0)goto 17071                                         2735
      eqs=max(eps,flmin)                                                   2735
      alf=eqs**(1.0/(nlam-1))                                              2735
17071 continue                                                             2736
      m=0                                                                  2736
      mm=0                                                                 2736
      nlp=0                                                                2736
      nin=nlp                                                              2736
      mnl=min(mnlam,nlam)                                                  2736
      as=0.0                                                               2736
      cthr=cthri*dev0                                                      2737
17080 do 17081 j=1,ni                                                      2737
      if(ju(j).eq.0)goto 17081                                             2737
      ga(j)=abs(dot_product(wr,x(:,j)))                                    2737
17081 continue                                                             2738
17082 continue                                                             2738
17090 do 17091 ilm=1,nlam                                                  2738
      al0=al                                                               2739
      if(flmin .lt. 1.0)goto 17111                                         2739
      al=ulam(ilm)                                                         2739
      goto 17101                                                           2740
17111 if(ilm .le. 2)goto 17121                                             2740
      al=al*alf                                                            2740
      goto 17101                                                           2741
17121 if(ilm .ne. 1)goto 17131                                             2741
      al=big                                                               2741
      goto 17141                                                           2742
17131 continue                                                             2742
      al0=0.0                                                              2743
17150 do 17151 j=1,ni                                                      2743
      if(ju(j).eq.0)goto 17151                                             2743
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2743
17151 continue                                                             2744
17152 continue                                                             2744
      al0=al0/max(parm,1.0d-3)                                             2744
      al=alf*al0                                                           2745
17141 continue                                                             2746
17101 continue                                                             2746
      sa=alpha*al                                                          2746
      omal=oma*al                                                          2746
      tlam=alpha*(2.0*al-al0)                                              2747
17160 do 17161 k=1,ni                                                      2747
      if(ixx(k).eq.1)goto 17161                                            2747
      if(ju(k).eq.0)goto 17161                                             2748
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     2749
17161 continue                                                             2750
17162 continue                                                             2750
11020 continue                                                             2751
17170 continue                                                             2751
17171 continue                                                             2751
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2752
      call vars(no,ni,x,w,ixx,v)                                           2753
17180 continue                                                             2753
17181 continue                                                             2753
      nlp=nlp+1                                                            2753
      dli=0.0                                                              2754
17190 do 17191 j=1,ni                                                      2754
      if(ixx(j).eq.0)goto 17191                                            2755
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                   2756
      if(abs(u) .gt. vp(j)*sa)goto 17211                                   2756
      at=0.0                                                               2756
      goto 17221                                                           2757
17211 continue                                                             2757
      at=max(cl(1,j),min(cl(2,j),sign(abs(u)-vp(j)*sa,u)/  (v(j)+vp(j)*o   2759 
     *mal)))
17221 continue                                                             2760
17201 continue                                                             2760
      if(at .eq. a(j))goto 17241                                           2760
      del=at-a(j)                                                          2760
      a(j)=at                                                              2760
      dli=max(dli,v(j)*del**2)                                             2761
      wr=wr-del*w*x(:,j)                                                   2761
      f=f+del*x(:,j)                                                       2762
      if(mm(j) .ne. 0)goto 17261                                           2762
      nin=nin+1                                                            2762
      if(nin.gt.nx)goto 17192                                              2763
      mm(j)=nin                                                            2763
      m(nin)=j                                                             2764
17261 continue                                                             2765
17241 continue                                                             2766
17191 continue                                                             2767
17192 continue                                                             2767
      if(nin.gt.nx)goto 17182                                              2767
      if(dli.lt.cthr)goto 17182                                            2768
      if(nlp .le. maxit)goto 17281                                         2768
      jerr=-ilm                                                            2768
      return                                                               2768
17281 continue                                                             2769
17290 continue                                                             2769
17291 continue                                                             2769
      nlp=nlp+1                                                            2769
      dli=0.0                                                              2770
17300 do 17301 l=1,nin                                                     2770
      j=m(l)                                                               2771
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                   2772
      if(abs(u) .gt. vp(j)*sa)goto 17321                                   2772
      at=0.0                                                               2772
      goto 17331                                                           2773
17321 continue                                                             2773
      at=max(cl(1,j),min(cl(2,j),sign(abs(u)-vp(j)*sa,u)/  (v(j)+vp(j)*o   2775 
     *mal)))
17331 continue                                                             2776
17311 continue                                                             2776
      if(at .eq. a(j))goto 17351                                           2776
      del=at-a(j)                                                          2776
      a(j)=at                                                              2776
      dli=max(dli,v(j)*del**2)                                             2777
      wr=wr-del*w*x(:,j)                                                   2777
      f=f+del*x(:,j)                                                       2778
17351 continue                                                             2779
17301 continue                                                             2780
17302 continue                                                             2780
      if(dli.lt.cthr)goto 17292                                            2780
      if(nlp .le. maxit)goto 17371                                         2780
      jerr=-ilm                                                            2780
      return                                                               2780
17371 continue                                                             2781
      goto 17291                                                           2782
17292 continue                                                             2782
      goto 17181                                                           2783
17182 continue                                                             2783
      if(nin.gt.nx)goto 17172                                              2784
      e=q*exp(sign(min(abs(f),fmax),f))                                    2785
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                         2786
      if(jerr .eq. 0)goto 17391                                            2786
      jerr=jerr-ilm                                                        2786
      go to 12320                                                          2786
17391 continue                                                             2787
      ix=0                                                                 2788
17400 do 17401 j=1,nin                                                     2788
      k=m(j)                                                               2789
      if(v(k)*(a(k)-as(k))**2.lt.cthr)goto 17401                           2789
      ix=1                                                                 2789
      goto 17402                                                           2789
17401 continue                                                             2790
17402 continue                                                             2790
      if(ix .ne. 0)goto 17421                                              2791
17430 do 17431 k=1,ni                                                      2791
      if(ixx(k).eq.1)goto 17431                                            2791
      if(ju(k).eq.0)goto 17431                                             2792
      ga(k)=abs(dot_product(wr,x(:,k)))                                    2793
      if(ga(k) .le. sa*vp(k))goto 17451                                    2793
      ixx(k)=1                                                             2793
      ix=1                                                                 2793
17451 continue                                                             2794
17431 continue                                                             2795
17432 continue                                                             2795
      if(ix.eq.1) go to 11020                                              2796
      goto 17172                                                           2797
17421 continue                                                             2798
      goto 17171                                                           2799
17172 continue                                                             2799
      if(nin .le. nx)goto 17471                                            2799
      jerr=-10000-ilm                                                      2799
      goto 17092                                                           2799
17471 continue                                                             2800
      if(nin.gt.0) ao(1:nin,ilm)=a(m(1:nin))                               2800
      kin(ilm)=nin                                                         2801
      alm(ilm)=al                                                          2801
      lmu=ilm                                                              2802
      dev(ilm)=(risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)-r0)/rr                   2803
      if(ilm.lt.mnl)goto 17091                                             2803
      if(flmin.ge.1.0)goto 17091                                           2804
      me=0                                                                 2804
17480 do 17481 j=1,nin                                                     2804
      if(ao(j,ilm).ne.0.0) me=me+1                                         2804
17481 continue                                                             2804
17482 continue                                                             2804
      if(me.gt.ne)goto 17092                                               2805
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 17092              2806
      if(dev(ilm).gt.devmax)goto 17092                                     2807
17091 continue                                                             2808
17092 continue                                                             2808
      g=f                                                                  2809
12320 continue                                                             2809
      deallocate(e,uu,w,dk,v,xs,f,wr,a,as,jp,kp,dq,mm,ga,ixx)              2810
      return                                                               2811
      end                                                                  2812
      subroutine cxmodval(ca,ia,nin,n,x,f)                                 2813
      implicit double precision(a-h,o-z)                                   2814
      double precision ca(nin),x(n,*),f(n)                                 2814
      integer ia(nin)                                                      2815
      f=0.0                                                                2815
      if(nin.le.0) return                                                  2816
17490 do 17491 i=1,n                                                       2816
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                      2816
17491 continue                                                             2817
17492 continue                                                             2817
      return                                                               2818
      end                                                                  2819
      subroutine groups(no,y,d,q,nk,kp,jp,t0,jerr)                         2820
      implicit double precision(a-h,o-z)                                   2821
      double precision y(no),d(no),q(no)                                   2821
      integer jp(no),kp(*)                                                 2822
17500 do 17501 j=1,no                                                      2822
      jp(j)=j                                                              2822
17501 continue                                                             2822
17502 continue                                                             2822
      call psort7(y,jp,1,no)                                               2823
      nj=0                                                                 2823
17510 do 17511 j=1,no                                                      2823
      if(q(jp(j)).le.0.0)goto 17511                                        2823
      nj=nj+1                                                              2823
      jp(nj)=jp(j)                                                         2823
17511 continue                                                             2824
17512 continue                                                             2824
      if(nj .ne. 0)goto 17531                                              2824
      jerr=20000                                                           2824
      return                                                               2824
17531 continue                                                             2825
      j=1                                                                  2825
17540 continue                                                             2825
17541 if(d(jp(j)).gt.0.0)goto 17542                                        2825
      j=j+1                                                                2825
      if(j.gt.nj)goto 17542                                                2825
      goto 17541                                                           2826
17542 continue                                                             2826
      if(j .lt. nj-1)goto 17561                                            2826
      jerr=30000                                                           2826
      return                                                               2826
17561 continue                                                             2827
      t0=y(jp(j))                                                          2827
      j0=j-1                                                               2828
      if(j0 .le. 0)goto 17581                                              2829
17590 continue                                                             2829
17591 if(y(jp(j0)).lt.t0)goto 17592                                        2829
      j0=j0-1                                                              2829
      if(j0.eq.0)goto 17592                                                2829
      goto 17591                                                           2830
17592 continue                                                             2830
      if(j0 .le. 0)goto 17611                                              2830
      nj=nj-j0                                                             2830
17620 do 17621 j=1,nj                                                      2830
      jp(j)=jp(j+j0)                                                       2830
17621 continue                                                             2830
17622 continue                                                             2830
17611 continue                                                             2831
17581 continue                                                             2832
      jerr=0                                                               2832
      nk=0                                                                 2832
      yk=t0                                                                2832
      j=2                                                                  2833
17630 continue                                                             2833
17631 continue                                                             2833
17640 continue                                                             2834
17641 if(d(jp(j)).gt.0.0.and.y(jp(j)).gt.yk)goto 17642                     2834
      j=j+1                                                                2834
      if(j.gt.nj)goto 17642                                                2834
      goto 17641                                                           2835
17642 continue                                                             2835
      nk=nk+1                                                              2835
      kp(nk)=j-1                                                           2835
      if(j.gt.nj)goto 17632                                                2836
      if(j .ne. nj)goto 17661                                              2836
      nk=nk+1                                                              2836
      kp(nk)=nj                                                            2836
      goto 17632                                                           2836
17661 continue                                                             2837
      yk=y(jp(j))                                                          2837
      j=j+1                                                                2838
      goto 17631                                                           2839
17632 continue                                                             2839
      return                                                               2840
      end                                                                  2841
      subroutine outer(no,nk,d,dk,kp,jp,e,wr,w,jerr,u)                     2842
      implicit double precision(a-h,o-z)                                   2843
      double precision d(no),dk(nk),wr(no),w(no)                           2844
      double precision e(no),u(no),b,c                                     2844
      integer kp(nk),jp(no)                                                2845
      call usk(no,nk,kp,jp,e,u)                                            2846
      b=dk(1)/u(1)                                                         2846
      c=dk(1)/u(1)**2                                                      2846
      jerr=0                                                               2847
17670 do 17671 j=1,kp(1)                                                   2847
      i=jp(j)                                                              2848
      w(i)=e(i)*(b-e(i)*c)                                                 2848
      if(w(i) .gt. 0.0)goto 17691                                          2848
      jerr=-30000                                                          2848
      return                                                               2848
17691 continue                                                             2849
      wr(i)=d(i)-e(i)*b                                                    2850
17671 continue                                                             2851
17672 continue                                                             2851
17700 do 17701 k=2,nk                                                      2851
      j1=kp(k-1)+1                                                         2851
      j2=kp(k)                                                             2852
      b=b+dk(k)/u(k)                                                       2852
      c=c+dk(k)/u(k)**2                                                    2853
17710 do 17711 j=j1,j2                                                     2853
      i=jp(j)                                                              2854
      w(i)=e(i)*(b-e(i)*c)                                                 2854
      if(w(i) .gt. 0.0)goto 17731                                          2854
      jerr=-30000                                                          2854
      return                                                               2854
17731 continue                                                             2855
      wr(i)=d(i)-e(i)*b                                                    2856
17711 continue                                                             2857
17712 continue                                                             2857
17701 continue                                                             2858
17702 continue                                                             2858
      return                                                               2859
      end                                                                  2860
      subroutine vars(no,ni,x,w,ixx,v)                                     2861
      implicit double precision(a-h,o-z)                                   2862
      double precision x(no,ni),w(no),v(ni)                                2862
      integer ixx(ni)                                                      2863
17740 do 17741 j=1,ni                                                      2863
      if(ixx(j).gt.0) v(j)=dot_product(w,x(:,j)**2)                        2863
17741 continue                                                             2864
17742 continue                                                             2864
      return                                                               2865
      end                                                                  2866
      subroutine died(no,nk,d,kp,jp,dk)                                    2867
      implicit double precision(a-h,o-z)                                   2868
      double precision d(no),dk(nk)                                        2868
      integer kp(nk),jp(no)                                                2869
      dk(1)=sum(d(jp(1:kp(1))))                                            2870
17750 do 17751 k=2,nk                                                      2870
      dk(k)=sum(d(jp((kp(k-1)+1):kp(k))))                                  2870
17751 continue                                                             2871
17752 continue                                                             2871
      return                                                               2872
      end                                                                  2873
      subroutine usk(no,nk,kp,jp,e,u)                                      2874
      implicit double precision(a-h,o-z)                                   2875
      double precision e(no),u(nk),h                                       2875
      integer kp(nk),jp(no)                                                2876
      h=0.0                                                                2877
17760 do 17761 k=nk,1,-1                                                   2877
      j2=kp(k)                                                             2878
      j1=1                                                                 2878
      if(k.gt.1) j1=kp(k-1)+1                                              2879
17770 do 17771 j=j2,j1,-1                                                  2879
      h=h+e(jp(j))                                                         2879
17771 continue                                                             2880
17772 continue                                                             2880
      u(k)=h                                                               2881
17761 continue                                                             2882
17762 continue                                                             2882
      return                                                               2883
      end                                                                  2884
      function risk(no,ni,nk,d,dk,f,e,kp,jp,u)                             2885
      implicit double precision(a-h,o-z)                                   2886
      double precision d(no),dk(nk),f(no)                                  2887
      integer kp(nk),jp(no)                                                2887
      double precision e(no),u(nk),s                                       2888
      call usk(no,nk,kp,jp,e,u)                                            2888
      u=log(u)                                                             2889
      risk=dot_product(d,f)-dot_product(dk,u)                              2890
      return                                                               2891
      end                                                                  2892
      subroutine loglike(no,ni,x,y,d,g,w,nlam,a,flog,jerr)                 2893
      implicit double precision(a-h,o-z)                                   2894
      double precision x(no,ni),y(no),d(no),g(no),w(no),a(ni,nlam),flog(   2895 
     *nlam)
      double precision, dimension (:), allocatable :: dk,f,xm,dq,q              
      double precision, dimension (:), allocatable :: e,uu                      
      integer, dimension (:), allocatable :: jp,kp                              
      allocate(e(1:no),stat=jerr)                                          2901
      if(jerr.ne.0) go to 12320                                            2902
      allocate(q(1:no),stat=jerr)                                          2903
      if(jerr.ne.0) go to 12320                                            2904
      allocate(uu(1:no),stat=jerr)                                         2905
      if(jerr.ne.0) go to 12320                                            2906
      allocate(f(1:no),stat=jerr)                                          2907
      if(jerr.ne.0) go to 12320                                            2908
      allocate(dk(1:no),stat=jerr)                                         2909
      if(jerr.ne.0) go to 12320                                            2910
      allocate(jp(1:no),stat=jerr)                                         2911
      if(jerr.ne.0) go to 12320                                            2912
      allocate(kp(1:no),stat=jerr)                                         2913
      if(jerr.ne.0) go to 12320                                            2914
      allocate(dq(1:no),stat=jerr)                                         2915
      if(jerr.ne.0) go to 12320                                            2916
      allocate(xm(1:ni),stat=jerr)                                         2917
      if(jerr.ne.0) go to 12320                                            2918
      q=max(0d0,w)                                                         2918
      sw=sum(q)                                                            2919
      if(sw .gt. 0.0)goto 17791                                            2919
      jerr=9999                                                            2919
      go to 12320                                                          2919
17791 continue                                                             2920
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                               2921
      if(jerr.ne.0) go to 12320                                            2921
      fmax=log(huge(e(1))*0.1)                                             2922
      dq=d*q                                                               2922
      call died(no,nk,dq,kp,jp,dk)                                         2922
      gm=dot_product(q,g)/sw                                               2923
17800 do 17801 j=1,ni                                                      2923
      xm(j)=dot_product(q,x(:,j))/sw                                       2923
17801 continue                                                             2924
17802 continue                                                             2924
17810 do 17811 lam=1,nlam                                                  2925
17820 do 17821 i=1,no                                                      2925
      f(i)=g(i)-gm+dot_product(a(:,lam),(x(i,:)-xm))                       2926
      e(i)=q(i)*exp(sign(min(abs(f(i)),fmax),f(i)))                        2927
17821 continue                                                             2928
17822 continue                                                             2928
      flog(lam)=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                          2929
17811 continue                                                             2930
17812 continue                                                             2930
12320 continue                                                             2930
      deallocate(e,uu,dk,f,jp,kp,dq)                                       2931
      return                                                               2932
      end                                                                  2933
      subroutine fishnet (parm,no,ni,x,y,g,w,jd,vp,cl,ne,nx,nlam,flmin,u   2935 
     *lam,thr,  isd,intr,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   2936
      double precision x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam)        2937
      double precision ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam),cl(2,ni)   2938
      integer jd(*),ia(nx),nin(nlam)                                       2939
      double precision, dimension (:), allocatable :: xm,xs,ww,vq               
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 17841                                    2943
      jerr=10000                                                           2943
      return                                                               2943
17841 continue                                                             2944
      if(minval(y) .ge. 0.0)goto 17861                                     2944
      jerr=8888                                                            2944
      return                                                               2944
17861 continue                                                             2945
      allocate(ww(1:no),stat=jerr)                                         2946
      if(jerr.ne.0) return                                                 2947
      allocate(ju(1:ni),stat=jerr)                                         2948
      if(jerr.ne.0) return                                                 2949
      allocate(vq(1:ni),stat=jerr)                                         2950
      if(jerr.ne.0) return                                                 2951
      allocate(xm(1:ni),stat=jerr)                                         2952
      if(jerr.ne.0) return                                                 2953
      if(isd .le. 0)goto 17881                                             2953
      allocate(xs(1:ni),stat=jerr)                                         2953
      if(jerr.ne.0) return                                                 2953
17881 continue                                                             2954
      call chkvars(no,ni,x,ju)                                             2955
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2956
      if(maxval(ju) .gt. 0)goto 17901                                      2956
      jerr=7777                                                            2956
      go to 12320                                                          2956
17901 continue                                                             2957
      vq=max(0d0,vp)                                                       2957
      vq=vq*ni/sum(vq)                                                     2958
      ww=max(0d0,w)                                                        2958
      sw=sum(ww)                                                           2958
      if(sw .gt. 0.0)goto 17921                                            2958
      jerr=9999                                                            2958
      go to 12320                                                          2958
17921 continue                                                             2959
      ww=ww/sw                                                             2960
      call lstandard1(no,ni,x,ww,ju,isd,intr,xm,xs)                        2961
      if(isd .le. 0)goto 17941                                             2961
17950 do 17951 j=1,ni                                                      2961
      cl(:,j)=cl(:,j)*xs(j)                                                2961
17951 continue                                                             2961
17952 continue                                                             2961
17941 continue                                                             2962
      call fishnet1(parm,no,ni,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam,t   2964 
     *hr,  isd,intr,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) go to 12320                                            2964
      dev0=2.0*sw*dev0                                                     2965
17960 do 17961 k=1,lmu                                                     2965
      nk=nin(k)                                                            2966
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                      2967
      if(intr .ne. 0)goto 17981                                            2967
      a0(k)=0.0                                                            2967
      goto 17991                                                           2968
17981 continue                                                             2968
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                     2968
17991 continue                                                             2969
17971 continue                                                             2969
17961 continue                                                             2970
17962 continue                                                             2970
12320 continue                                                             2970
      deallocate(ww,ju,vq,xm)                                              2970
      if(isd.gt.0) deallocate(xs)                                          2971
      return                                                               2972
      end                                                                  2973
      subroutine fishnet1(parm,no,ni,x,y,g,q,ju,vp,cl,ne,nx,nlam,flmin,u   2975 
     *lam,shri,  isd,intr,maxit,lmu,a0,ca,m,kin,dev0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   2976
      double precision x(no,ni),y(no),g(no),q(no),vp(ni),ulam(nlam)        2977
      double precision ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam),cl(2,ni)   2978
      integer ju(ni),m(nx),kin(nlam)                                       2979
      double precision, dimension (:), allocatable :: t,w,wr,v,a,f,as,ga        
      integer, dimension (:), allocatable :: mm,ixx                             
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               2983
      sml=sml*10.0                                                         2984
      allocate(a(1:ni),stat=jerr)                                          2985
      if(jerr.ne.0) return                                                 2986
      allocate(as(1:ni),stat=jerr)                                         2987
      if(jerr.ne.0) return                                                 2988
      allocate(t(1:no),stat=jerr)                                          2989
      if(jerr.ne.0) return                                                 2990
      allocate(mm(1:ni),stat=jerr)                                         2991
      if(jerr.ne.0) return                                                 2992
      allocate(ga(1:ni),stat=jerr)                                         2993
      if(jerr.ne.0) return                                                 2994
      allocate(ixx(1:ni),stat=jerr)                                        2995
      if(jerr.ne.0) return                                                 2996
      allocate(wr(1:no),stat=jerr)                                         2997
      if(jerr.ne.0) return                                                 2998
      allocate(v(1:ni),stat=jerr)                                          2999
      if(jerr.ne.0) return                                                 3000
      allocate(w(1:no),stat=jerr)                                          3001
      if(jerr.ne.0) return                                                 3002
      allocate(f(1:no),stat=jerr)                                          3003
      if(jerr.ne.0) return                                                 3004
      bta=parm                                                             3004
      omb=1.0-bta                                                          3005
      t=q*y                                                                3005
      yb=sum(t)                                                            3005
      fmax=log(huge(bta)*0.1)                                              3006
      if(nonzero(no,g) .ne. 0)goto 18011                                   3007
      if(intr .eq. 0)goto 18031                                            3007
      w=q*yb                                                               3007
      az=log(yb)                                                           3007
      f=az                                                                 3007
      dv0=yb*(az-1.0)                                                      3007
      goto 18041                                                           3008
18031 continue                                                             3008
      w=q                                                                  3008
      az=0.0                                                               3008
      f=az                                                                 3008
      dv0=-1.0                                                             3008
18041 continue                                                             3009
18021 continue                                                             3009
      goto 18051                                                           3010
18011 continue                                                             3010
      w=q*exp(sign(min(abs(g),fmax),g))                                    3010
      v0=sum(w)                                                            3011
      if(intr .eq. 0)goto 18071                                            3011
      eaz=yb/v0                                                            3011
      w=eaz*w                                                              3011
      az=log(eaz)                                                          3011
      f=az+g                                                               3012
      dv0=dot_product(t,g)-yb*(1.0-az)                                     3013
      goto 18081                                                           3014
18071 continue                                                             3014
      az=0.0                                                               3014
      f=g                                                                  3014
      dv0=dot_product(t,g)-v0                                              3014
18081 continue                                                             3015
18061 continue                                                             3015
18051 continue                                                             3016
18001 continue                                                             3016
      a=0.0                                                                3016
      as=0.0                                                               3016
      wr=t-w                                                               3016
      v0=1.0                                                               3016
      if(intr.ne.0) v0=yb                                                  3016
      dvr=-yb                                                              3017
18090 do 18091 i=1,no                                                      3017
      if(t(i).gt.0.0) dvr=dvr+t(i)*log(y(i))                               3017
18091 continue                                                             3017
18092 continue                                                             3017
      dvr=dvr-dv0                                                          3017
      dev0=dvr                                                             3019
      alf=1.0                                                              3021
      if(flmin .ge. 1.0)goto 18111                                         3021
      eqs=max(eps,flmin)                                                   3021
      alf=eqs**(1.0/(nlam-1))                                              3021
18111 continue                                                             3022
      m=0                                                                  3022
      mm=0                                                                 3022
      nlp=0                                                                3022
      nin=nlp                                                              3022
      mnl=min(mnlam,nlam)                                                  3022
      shr=shri*dev0                                                        3022
      ixx=0                                                                3022
      al=0.0                                                               3023
18120 do 18121 j=1,ni                                                      3023
      if(ju(j).eq.0)goto 18121                                             3023
      ga(j)=abs(dot_product(wr,x(:,j)))                                    3023
18121 continue                                                             3024
18122 continue                                                             3024
18130 do 18131 ilm=1,nlam                                                  3024
      al0=al                                                               3025
      if(flmin .lt. 1.0)goto 18151                                         3025
      al=ulam(ilm)                                                         3025
      goto 18141                                                           3026
18151 if(ilm .le. 2)goto 18161                                             3026
      al=al*alf                                                            3026
      goto 18141                                                           3027
18161 if(ilm .ne. 1)goto 18171                                             3027
      al=big                                                               3027
      goto 18181                                                           3028
18171 continue                                                             3028
      al0=0.0                                                              3029
18190 do 18191 j=1,ni                                                      3029
      if(ju(j).eq.0)goto 18191                                             3029
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            3029
18191 continue                                                             3030
18192 continue                                                             3030
      al0=al0/max(bta,1.0d-3)                                              3030
      al=alf*al0                                                           3031
18181 continue                                                             3032
18141 continue                                                             3032
      al2=al*omb                                                           3032
      al1=al*bta                                                           3032
      tlam=bta*(2.0*al-al0)                                                3033
18200 do 18201 k=1,ni                                                      3033
      if(ixx(k).eq.1)goto 18201                                            3033
      if(ju(k).eq.0)goto 18201                                             3034
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     3035
18201 continue                                                             3036
18202 continue                                                             3036
11020 continue                                                             3037
18210 continue                                                             3037
18211 continue                                                             3037
      az0=az                                                               3038
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                3039
18220 do 18221 j=1,ni                                                      3039
      if(ixx(j).ne.0) v(j)=dot_product(w,x(:,j)**2)                        3039
18221 continue                                                             3040
18222 continue                                                             3040
18230 continue                                                             3040
18231 continue                                                             3040
      nlp=nlp+1                                                            3040
      dlx=0.0                                                              3041
18240 do 18241 k=1,ni                                                      3041
      if(ixx(k).eq.0)goto 18241                                            3041
      ak=a(k)                                                              3042
      u=dot_product(wr,x(:,k))+v(k)*ak                                     3042
      au=abs(u)-vp(k)*al1                                                  3043
      if(au .gt. 0.0)goto 18261                                            3043
      a(k)=0.0                                                             3043
      goto 18271                                                           3044
18261 continue                                                             3045
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))           3046
18271 continue                                                             3047
18251 continue                                                             3047
      if(a(k).eq.ak)goto 18241                                             3047
      d=a(k)-ak                                                            3047
      dlx=max(dlx,v(k)*d**2)                                               3048
      wr=wr-d*w*x(:,k)                                                     3048
      f=f+d*x(:,k)                                                         3049
      if(mm(k) .ne. 0)goto 18291                                           3049
      nin=nin+1                                                            3049
      if(nin.gt.nx)goto 18242                                              3050
      mm(k)=nin                                                            3050
      m(nin)=k                                                             3051
18291 continue                                                             3052
18241 continue                                                             3053
18242 continue                                                             3053
      if(nin.gt.nx)goto 18232                                              3054
      if(intr .eq. 0)goto 18311                                            3054
      d=sum(wr)/v0                                                         3055
      az=az+d                                                              3055
      dlx=max(dlx,v0*d**2)                                                 3055
      wr=wr-d*w                                                            3055
      f=f+d                                                                3056
18311 continue                                                             3057
      if(dlx.lt.shr)goto 18232                                             3057
      if(nlp .le. maxit)goto 18331                                         3057
      jerr=-ilm                                                            3057
      return                                                               3057
18331 continue                                                             3058
18340 continue                                                             3058
18341 continue                                                             3058
      nlp=nlp+1                                                            3058
      dlx=0.0                                                              3059
18350 do 18351 l=1,nin                                                     3059
      k=m(l)                                                               3059
      ak=a(k)                                                              3060
      u=dot_product(wr,x(:,k))+v(k)*ak                                     3060
      au=abs(u)-vp(k)*al1                                                  3061
      if(au .gt. 0.0)goto 18371                                            3061
      a(k)=0.0                                                             3061
      goto 18381                                                           3062
18371 continue                                                             3063
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))           3064
18381 continue                                                             3065
18361 continue                                                             3065
      if(a(k).eq.ak)goto 18351                                             3065
      d=a(k)-ak                                                            3065
      dlx=max(dlx,v(k)*d**2)                                               3066
      wr=wr-d*w*x(:,k)                                                     3066
      f=f+d*x(:,k)                                                         3068
18351 continue                                                             3068
18352 continue                                                             3068
      if(intr .eq. 0)goto 18401                                            3068
      d=sum(wr)/v0                                                         3068
      az=az+d                                                              3069
      dlx=max(dlx,v0*d**2)                                                 3069
      wr=wr-d*w                                                            3069
      f=f+d                                                                3070
18401 continue                                                             3071
      if(dlx.lt.shr)goto 18342                                             3071
      if(nlp .le. maxit)goto 18421                                         3071
      jerr=-ilm                                                            3071
      return                                                               3071
18421 continue                                                             3072
      goto 18341                                                           3073
18342 continue                                                             3073
      goto 18231                                                           3074
18232 continue                                                             3074
      if(nin.gt.nx)goto 18212                                              3075
      w=q*exp(sign(min(abs(f),fmax),f))                                    3075
      v0=sum(w)                                                            3075
      wr=t-w                                                               3076
      if(v0*(az-az0)**2 .ge. shr)goto 18441                                3076
      ix=0                                                                 3077
18450 do 18451 j=1,nin                                                     3077
      k=m(j)                                                               3078
      if(v(k)*(a(k)-as(k))**2.lt.shr)goto 18451                            3078
      ix=1                                                                 3078
      goto 18452                                                           3079
18451 continue                                                             3080
18452 continue                                                             3080
      if(ix .ne. 0)goto 18471                                              3081
18480 do 18481 k=1,ni                                                      3081
      if(ixx(k).eq.1)goto 18481                                            3081
      if(ju(k).eq.0)goto 18481                                             3082
      ga(k)=abs(dot_product(wr,x(:,k)))                                    3083
      if(ga(k) .le. al1*vp(k))goto 18501                                   3083
      ixx(k)=1                                                             3083
      ix=1                                                                 3083
18501 continue                                                             3084
18481 continue                                                             3085
18482 continue                                                             3085
      if(ix.eq.1) go to 11020                                              3086
      goto 18212                                                           3087
18471 continue                                                             3088
18441 continue                                                             3089
      goto 18211                                                           3090
18212 continue                                                             3090
      if(nin .le. nx)goto 18521                                            3090
      jerr=-10000-ilm                                                      3090
      goto 18132                                                           3090
18521 continue                                                             3091
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                               3091
      kin(ilm)=nin                                                         3092
      a0(ilm)=az                                                           3092
      alm(ilm)=al                                                          3092
      lmu=ilm                                                              3093
      dev(ilm)=(dot_product(t,f)-v0-dv0)/dvr                               3094
      if(ilm.lt.mnl)goto 18131                                             3094
      if(flmin.ge.1.0)goto 18131                                           3095
      me=0                                                                 3095
18530 do 18531 j=1,nin                                                     3095
      if(ca(j,ilm).ne.0.0) me=me+1                                         3095
18531 continue                                                             3095
18532 continue                                                             3095
      if(me.gt.ne)goto 18132                                               3096
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 18132              3097
      if(dev(ilm).gt.devmax)goto 18132                                     3098
18131 continue                                                             3099
18132 continue                                                             3099
      g=f                                                                  3100
12320 continue                                                             3100
      deallocate(t,w,wr,v,a,f,as,mm,ga,ixx)                                3101
      return                                                               3102
      end                                                                  3103
      function nonzero(n,v)                                                3104
      implicit double precision(a-h,o-z)                                   3105
      double precision v(n)                                                3106
      nonzero=0                                                            3106
18540 do 18541 i=1,n                                                       3106
      if(v(i) .eq. 0.0)goto 18561                                          3106
      nonzero=1                                                            3106
      return                                                               3106
18561 continue                                                             3106
18541 continue                                                             3107
18542 continue                                                             3107
      return                                                               3108
      end                                                                  3109
      subroutine solns(ni,nx,lmu,a,ia,nin,b)                               3110
      implicit double precision(a-h,o-z)                                   3111
      double precision a(nx,lmu),b(ni,lmu)                                 3111
      integer ia(nx),nin(lmu)                                              3112
18570 do 18571 lam=1,lmu                                                   3112
      call uncomp(ni,a(:,lam),ia,nin(lam),b(:,lam))                        3112
18571 continue                                                             3113
18572 continue                                                             3113
      return                                                               3114
      end                                                                  3115
      subroutine lsolns(ni,nx,nc,lmu,a,ia,nin,b)                           3116
      implicit double precision(a-h,o-z)                                   3117
      double precision a(nx,nc,lmu),b(ni,nc,lmu)                           3117
      integer ia(nx),nin(lmu)                                              3118
18580 do 18581 lam=1,lmu                                                   3118
      call luncomp(ni,nx,nc,a(1,1,lam),ia,nin(lam),b(1,1,lam))             3118
18581 continue                                                             3119
18582 continue                                                             3119
      return                                                               3120
      end                                                                  3121
      subroutine deviance(no,ni,x,y,g,q,nlam,a0,a,flog,jerr)               3122
      implicit double precision(a-h,o-z)                                   3123
      double precision x(no,ni),y(no),g(no),q(no),a(ni,nlam),a0(nlam),fl   3124 
     *og(nlam)
      double precision, dimension (:), allocatable :: w                         
      if(minval(y) .ge. 0.0)goto 18601                                     3127
      jerr=8888                                                            3127
      return                                                               3127
18601 continue                                                             3128
      allocate(w(1:no),stat=jerr)                                          3128
      if(jerr.ne.0) return                                                 3129
      w=max(0d0,q)                                                         3129
      sw=sum(w)                                                            3129
      if(sw .gt. 0.0)goto 18621                                            3129
      jerr=9999                                                            3129
      go to 12320                                                          3129
18621 continue                                                             3130
      yb=dot_product(w,y)/sw                                               3130
      fmax=log(huge(y(1))*0.1)                                             3131
18630 do 18631 lam=1,nlam                                                  3131
      s=0.0                                                                3132
18640 do 18641 i=1,no                                                      3132
      if(w(i).le.0.0)goto 18641                                            3133
      f=g(i)+a0(lam)+dot_product(a(:,lam),x(i,:))                          3134
      s=s+w(i)*(y(i)*f-exp(sign(min(abs(f),fmax),f)))                      3135
18641 continue                                                             3136
18642 continue                                                             3136
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                3137
18631 continue                                                             3138
18632 continue                                                             3138
12320 continue                                                             3138
      deallocate(w)                                                        3139
      return                                                               3140
      end                                                                  3141
      subroutine spfishnet (parm,no,ni,x,ix,jx,y,g,w,jd,vp,cl,ne,nx,nlam   3143 
     *,flmin,  ulam,thr,isd,intr,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp
     *,jerr)
      implicit double precision(a-h,o-z)                                   3144
      double precision x(*),y(no),g(no),w(no),vp(ni),ulam(nlam),cl(2,ni)   3145
      double precision ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)            3146
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           3147
      double precision, dimension (:), allocatable :: xm,xs,ww,vq               
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 18661                                    3151
      jerr=10000                                                           3151
      return                                                               3151
18661 continue                                                             3152
      if(minval(y) .ge. 0.0)goto 18681                                     3152
      jerr=8888                                                            3152
      return                                                               3152
18681 continue                                                             3153
      allocate(ww(1:no),stat=jerr)                                         3154
      if(jerr.ne.0) return                                                 3155
      allocate(ju(1:ni),stat=jerr)                                         3156
      if(jerr.ne.0) return                                                 3157
      allocate(vq(1:ni),stat=jerr)                                         3158
      if(jerr.ne.0) return                                                 3159
      allocate(xm(1:ni),stat=jerr)                                         3160
      if(jerr.ne.0) return                                                 3161
      allocate(xs(1:ni),stat=jerr)                                         3162
      if(jerr.ne.0) return                                                 3163
      call spchkvars(no,ni,x,ix,ju)                                        3164
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 3165
      if(maxval(ju) .gt. 0)goto 18701                                      3165
      jerr=7777                                                            3165
      go to 12320                                                          3165
18701 continue                                                             3166
      vq=max(0d0,vp)                                                       3166
      vq=vq*ni/sum(vq)                                                     3167
      ww=max(0d0,w)                                                        3167
      sw=sum(ww)                                                           3167
      if(sw .gt. 0.0)goto 18721                                            3167
      jerr=9999                                                            3167
      go to 12320                                                          3167
18721 continue                                                             3168
      ww=ww/sw                                                             3169
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs)                3170
      if(isd .le. 0)goto 18741                                             3170
18750 do 18751 j=1,ni                                                      3170
      cl(:,j)=cl(:,j)*xs(j)                                                3170
18751 continue                                                             3170
18752 continue                                                             3170
18741 continue                                                             3171
      call spfishnet1(parm,no,ni,x,ix,jx,y,g,ww,ju,vq,cl,ne,nx,nlam,flmi   3173 
     *n,ulam,thr,  isd,intr,maxit,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm,nl
     *p,jerr)
      if(jerr.gt.0) go to 12320                                            3173
      dev0=2.0*sw*dev0                                                     3174
18760 do 18761 k=1,lmu                                                     3174
      nk=nin(k)                                                            3175
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                      3176
      if(intr .ne. 0)goto 18781                                            3176
      a0(k)=0.0                                                            3176
      goto 18791                                                           3177
18781 continue                                                             3177
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                     3177
18791 continue                                                             3178
18771 continue                                                             3178
18761 continue                                                             3179
18762 continue                                                             3179
12320 continue                                                             3179
      deallocate(ww,ju,vq,xm,xs)                                           3180
      return                                                               3181
      end                                                                  3182
      subroutine spfishnet1(parm,no,ni,x,ix,jx,y,g,q,ju,vp,cl,ne,nx,nlam   3184 
     *,flmin,ulam,  shri,isd,intr,maxit,xb,xs,lmu,a0,ca,m,kin,dev0,dev,a
     *lm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   3185
      double precision x(*),y(no),g(no),q(no),vp(ni),ulam(nlam),xb(ni),x   3186 
     *s(ni)
      double precision ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam),cl(2,ni)   3187
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           3188
      double precision, dimension (:), allocatable :: qy,t,w,wr,v               
      double precision, dimension (:), allocatable :: a,as,xm,ga                
      integer, dimension (:), allocatable :: mm,ixx                             
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               3193
      sml=sml*10.0                                                         3194
      allocate(a(1:ni),stat=jerr)                                          3195
      if(jerr.ne.0) return                                                 3196
      allocate(as(1:ni),stat=jerr)                                         3197
      if(jerr.ne.0) return                                                 3198
      allocate(t(1:no),stat=jerr)                                          3199
      if(jerr.ne.0) return                                                 3200
      allocate(mm(1:ni),stat=jerr)                                         3201
      if(jerr.ne.0) return                                                 3202
      allocate(ga(1:ni),stat=jerr)                                         3203
      if(jerr.ne.0) return                                                 3204
      allocate(ixx(1:ni),stat=jerr)                                        3205
      if(jerr.ne.0) return                                                 3206
      allocate(wr(1:no),stat=jerr)                                         3207
      if(jerr.ne.0) return                                                 3208
      allocate(v(1:ni),stat=jerr)                                          3209
      if(jerr.ne.0) return                                                 3210
      allocate(xm(1:ni),stat=jerr)                                         3211
      if(jerr.ne.0) return                                                 3212
      allocate(w(1:no),stat=jerr)                                          3213
      if(jerr.ne.0) return                                                 3214
      allocate(qy(1:no),stat=jerr)                                         3215
      if(jerr.ne.0) return                                                 3216
      bta=parm                                                             3216
      omb=1.0-bta                                                          3216
      fmax=log(huge(bta)*0.1)                                              3217
      qy=q*y                                                               3217
      yb=sum(qy)                                                           3218
      if(nonzero(no,g) .ne. 0)goto 18811                                   3218
      t=0.0                                                                3219
      if(intr .eq. 0)goto 18831                                            3219
      w=q*yb                                                               3219
      az=log(yb)                                                           3219
      uu=az                                                                3220
      xm=yb*xb                                                             3220
      dv0=yb*(az-1.0)                                                      3221
      goto 18841                                                           3222
18831 continue                                                             3222
      w=q                                                                  3222
      xm=0.0                                                               3222
      uu=0.0                                                               3222
      az=uu                                                                3222
      dv0=-1.0                                                             3222
18841 continue                                                             3223
18821 continue                                                             3223
      goto 18851                                                           3224
18811 continue                                                             3224
      w=q*exp(sign(min(abs(g),fmax),g))                                    3224
      ww=sum(w)                                                            3224
      t=g                                                                  3225
      if(intr .eq. 0)goto 18871                                            3225
      eaz=yb/ww                                                            3226
      w=eaz*w                                                              3226
      az=log(eaz)                                                          3226
      uu=az                                                                3226
      dv0=dot_product(qy,g)-yb*(1.0-az)                                    3227
      goto 18881                                                           3228
18871 continue                                                             3228
      uu=0.0                                                               3228
      az=uu                                                                3228
      dv0=dot_product(qy,g)-ww                                             3228
18881 continue                                                             3229
18861 continue                                                             3229
18890 do 18891 j=1,ni                                                      3229
      if(ju(j).eq.0)goto 18891                                             3229
      jb=ix(j)                                                             3229
      je=ix(j+1)-1                                                         3230
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             3231
18891 continue                                                             3232
18892 continue                                                             3232
18851 continue                                                             3233
18801 continue                                                             3233
      tt=yb*uu                                                             3233
      ww=1.0                                                               3233
      if(intr.ne.0) ww=yb                                                  3233
      wr=qy-q*(yb*(1.0-uu))                                                3233
      a=0.0                                                                3233
      as=0.0                                                               3234
      dvr=-yb                                                              3235
18900 do 18901 i=1,no                                                      3235
      if(qy(i).gt.0.0) dvr=dvr+qy(i)*log(y(i))                             3235
18901 continue                                                             3235
18902 continue                                                             3235
      dvr=dvr-dv0                                                          3235
      dev0=dvr                                                             3237
      alf=1.0                                                              3239
      if(flmin .ge. 1.0)goto 18921                                         3239
      eqs=max(eps,flmin)                                                   3239
      alf=eqs**(1.0/(nlam-1))                                              3239
18921 continue                                                             3240
      m=0                                                                  3240
      mm=0                                                                 3240
      nlp=0                                                                3240
      nin=nlp                                                              3240
      mnl=min(mnlam,nlam)                                                  3240
      shr=shri*dev0                                                        3240
      al=0.0                                                               3240
      ixx=0                                                                3241
18930 do 18931 j=1,ni                                                      3241
      if(ju(j).eq.0)goto 18931                                             3242
      jb=ix(j)                                                             3242
      je=ix(j+1)-1                                                         3243
      ga(j)=abs(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(j)-ww*xb(j)   3245 
     *)-xb(j)*tt)/xs(j)
18931 continue                                                             3246
18932 continue                                                             3246
18940 do 18941 ilm=1,nlam                                                  3246
      al0=al                                                               3247
      if(flmin .lt. 1.0)goto 18961                                         3247
      al=ulam(ilm)                                                         3247
      goto 18951                                                           3248
18961 if(ilm .le. 2)goto 18971                                             3248
      al=al*alf                                                            3248
      goto 18951                                                           3249
18971 if(ilm .ne. 1)goto 18981                                             3249
      al=big                                                               3249
      goto 18991                                                           3250
18981 continue                                                             3250
      al0=0.0                                                              3251
19000 do 19001 j=1,ni                                                      3251
      if(ju(j).eq.0)goto 19001                                             3251
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            3251
19001 continue                                                             3252
19002 continue                                                             3252
      al0=al0/max(bta,1.0d-3)                                              3252
      al=alf*al0                                                           3253
18991 continue                                                             3254
18951 continue                                                             3254
      al2=al*omb                                                           3254
      al1=al*bta                                                           3254
      tlam=bta*(2.0*al-al0)                                                3255
19010 do 19011 k=1,ni                                                      3255
      if(ixx(k).eq.1)goto 19011                                            3255
      if(ju(k).eq.0)goto 19011                                             3256
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     3257
19011 continue                                                             3258
19012 continue                                                             3258
11020 continue                                                             3259
19020 continue                                                             3259
19021 continue                                                             3259
      az0=az                                                               3260
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                3261
19030 do 19031 j=1,ni                                                      3261
      if(ixx(j).eq.0)goto 19031                                            3261
      jb=ix(j)                                                             3261
      je=ix(j+1)-1                                                         3262
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             3263
      v(j)=(dot_product(w(jx(jb:je)),x(jb:je)**2)  -2.0*xb(j)*xm(j)+ww*x   3265 
     *b(j)**2)/xs(j)**2
19031 continue                                                             3266
19032 continue                                                             3266
19040 continue                                                             3266
19041 continue                                                             3266
      nlp=nlp+1                                                            3267
      dlx=0.0                                                              3268
19050 do 19051 k=1,ni                                                      3268
      if(ixx(k).eq.0)goto 19051                                            3268
      jb=ix(k)                                                             3268
      je=ix(k+1)-1                                                         3268
      ak=a(k)                                                              3269
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)   3271 
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                                  3272
      if(au .gt. 0.0)goto 19071                                            3272
      a(k)=0.0                                                             3272
      goto 19081                                                           3273
19071 continue                                                             3274
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))           3275
19081 continue                                                             3276
19061 continue                                                             3276
      if(a(k).eq.ak)goto 19051                                             3277
      if(mm(k) .ne. 0)goto 19101                                           3277
      nin=nin+1                                                            3277
      if(nin.gt.nx)goto 19052                                              3278
      mm(k)=nin                                                            3278
      m(nin)=k                                                             3279
19101 continue                                                             3280
      d=a(k)-ak                                                            3280
      dlx=max(dlx,v(k)*d**2)                                               3280
      dv=d/xs(k)                                                           3281
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)                 3282
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                                3283
      uu=uu-dv*xb(k)                                                       3283
      tt=tt-dv*xm(k)                                                       3284
19051 continue                                                             3285
19052 continue                                                             3285
      if(nin.gt.nx)goto 19042                                              3286
      if(intr .eq. 0)goto 19121                                            3286
      d=tt/ww-uu                                                           3287
      az=az+d                                                              3287
      dlx=max(dlx,ww*d**2)                                                 3287
      uu=uu+d                                                              3288
19121 continue                                                             3289
      if(dlx.lt.shr)goto 19042                                             3289
      if(nlp .le. maxit)goto 19141                                         3289
      jerr=-ilm                                                            3289
      return                                                               3289
19141 continue                                                             3290
19150 continue                                                             3290
19151 continue                                                             3290
      nlp=nlp+1                                                            3290
      dlx=0.0                                                              3291
19160 do 19161 l=1,nin                                                     3291
      k=m(l)                                                               3292
      jb=ix(k)                                                             3292
      je=ix(k+1)-1                                                         3292
      ak=a(k)                                                              3293
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)   3295 
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                                  3296
      if(au .gt. 0.0)goto 19181                                            3296
      a(k)=0.0                                                             3296
      goto 19191                                                           3297
19181 continue                                                             3298
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))           3299
19191 continue                                                             3300
19171 continue                                                             3300
      if(a(k).eq.ak)goto 19161                                             3300
      d=a(k)-ak                                                            3300
      dlx=max(dlx,v(k)*d**2)                                               3301
      dv=d/xs(k)                                                           3301
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)                 3302
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                                3303
      uu=uu-dv*xb(k)                                                       3303
      tt=tt-dv*xm(k)                                                       3304
19161 continue                                                             3305
19162 continue                                                             3305
      if(intr .eq. 0)goto 19211                                            3305
      d=tt/ww-uu                                                           3305
      az=az+d                                                              3306
      dlx=max(dlx,ww*d**2)                                                 3306
      uu=uu+d                                                              3307
19211 continue                                                             3308
      if(dlx.lt.shr)goto 19152                                             3308
      if(nlp .le. maxit)goto 19231                                         3308
      jerr=-ilm                                                            3308
      return                                                               3308
19231 continue                                                             3309
      goto 19151                                                           3310
19152 continue                                                             3310
      goto 19041                                                           3311
19042 continue                                                             3311
      if(nin.gt.nx)goto 19022                                              3312
      euu=exp(sign(min(abs(uu),fmax),uu))                                  3313
      w=euu*q*exp(sign(min(abs(t),fmax),t))                                3313
      ww=sum(w)                                                            3314
      wr=qy-w*(1.0-uu)                                                     3314
      tt=sum(wr)                                                           3315
      if(ww*(az-az0)**2 .ge. shr)goto 19251                                3315
      kx=0                                                                 3316
19260 do 19261 j=1,nin                                                     3316
      k=m(j)                                                               3317
      if(v(k)*(a(k)-as(k))**2.lt.shr)goto 19261                            3317
      kx=1                                                                 3317
      goto 19262                                                           3318
19261 continue                                                             3319
19262 continue                                                             3319
      if(kx .ne. 0)goto 19281                                              3320
19290 do 19291 j=1,ni                                                      3320
      if(ixx(j).eq.1)goto 19291                                            3320
      if(ju(j).eq.0)goto 19291                                             3321
      jb=ix(j)                                                             3321
      je=ix(j+1)-1                                                         3322
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             3323
      ga(j)=abs(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(j)-ww*xb(j)   3325 
     *)-xb(j)*tt)/xs(j)
      if(ga(j) .le. al1*vp(j))goto 19311                                   3325
      ixx(j)=1                                                             3325
      kx=1                                                                 3325
19311 continue                                                             3326
19291 continue                                                             3327
19292 continue                                                             3327
      if(kx.eq.1) go to 11020                                              3328
      goto 19022                                                           3329
19281 continue                                                             3330
19251 continue                                                             3331
      goto 19021                                                           3332
19022 continue                                                             3332
      if(nin .le. nx)goto 19331                                            3332
      jerr=-10000-ilm                                                      3332
      goto 18942                                                           3332
19331 continue                                                             3333
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                               3333
      kin(ilm)=nin                                                         3334
      a0(ilm)=az                                                           3334
      alm(ilm)=al                                                          3334
      lmu=ilm                                                              3335
      dev(ilm)=(dot_product(qy,t)+yb*uu-ww-dv0)/dvr                        3336
      if(ilm.lt.mnl)goto 18941                                             3336
      if(flmin.ge.1.0)goto 18941                                           3337
      me=0                                                                 3337
19340 do 19341 j=1,nin                                                     3337
      if(ca(j,ilm).ne.0.0) me=me+1                                         3337
19341 continue                                                             3337
19342 continue                                                             3337
      if(me.gt.ne)goto 18942                                               3338
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 18942              3339
      if(dev(ilm).gt.devmax)goto 18942                                     3340
18941 continue                                                             3341
18942 continue                                                             3341
      g=t+uu                                                               3342
12320 continue                                                             3342
      deallocate(t,w,wr,v,a,qy,xm,as,mm,ga,ixx)                            3343
      return                                                               3344
      end                                                                  3345
      subroutine spdeviance(no,ni,x,ix,jx,y,g,q,nlam,a0,a,flog,jerr)       3346
      implicit double precision(a-h,o-z)                                   3347
      double precision x(*),y(no),g(no),q(no),a(ni,nlam),a0(nlam),flog(n   3348 
     *lam)
      integer ix(*),jx(*)                                                  3349
      double precision, dimension (:), allocatable :: w,f                       
      if(minval(y) .ge. 0.0)goto 19361                                     3352
      jerr=8888                                                            3352
      return                                                               3352
19361 continue                                                             3353
      allocate(w(1:no),stat=jerr)                                          3354
      if(jerr.ne.0) return                                                 3355
      allocate(f(1:no),stat=jerr)                                          3356
      if(jerr.ne.0) return                                                 3357
      w=max(0d0,q)                                                         3357
      sw=sum(w)                                                            3357
      if(sw .gt. 0.0)goto 19381                                            3357
      jerr=9999                                                            3357
      go to 12320                                                          3357
19381 continue                                                             3358
      yb=dot_product(w,y)/sw                                               3358
      fmax=log(huge(y(1))*0.1)                                             3359
19390 do 19391 lam=1,nlam                                                  3359
      f=a0(lam)                                                            3360
19400 do 19401 j=1,ni                                                      3360
      if(a(j,lam).eq.0.0)goto 19401                                        3360
      jb=ix(j)                                                             3360
      je=ix(j+1)-1                                                         3361
      f(jx(jb:je))=f(jx(jb:je))+a(j,lam)*x(jb:je)                          3362
19401 continue                                                             3363
19402 continue                                                             3363
      f=f+g                                                                3364
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                   3365
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                3366
19391 continue                                                             3367
19392 continue                                                             3367
12320 continue                                                             3367
      deallocate(w,f)                                                      3368
      return                                                               3369
      end                                                                  3370
      subroutine cspdeviance(no,x,ix,jx,y,g,q,nx,nlam,a0,ca,ia,nin,flog,   3371 
     *jerr)
      implicit double precision(a-h,o-z)                                   3372
      double precision x(*),y(no),g(no),q(no),ca(nx,nlam),a0(nlam),flog(   3373 
     *nlam)
      integer ix(*),jx(*),nin(nlam),ia(nx)                                 3374
      double precision, dimension (:), allocatable :: w,f                       
      if(minval(y) .ge. 0.0)goto 19421                                     3377
      jerr=8888                                                            3377
      return                                                               3377
19421 continue                                                             3378
      allocate(w(1:no),stat=jerr)                                          3379
      if(jerr.ne.0) return                                                 3380
      allocate(f(1:no),stat=jerr)                                          3381
      if(jerr.ne.0) return                                                 3382
      w=max(0d0,q)                                                         3382
      sw=sum(w)                                                            3382
      if(sw .gt. 0.0)goto 19441                                            3382
      jerr=9999                                                            3382
      go to 12320                                                          3382
19441 continue                                                             3383
      yb=dot_product(w,y)/sw                                               3383
      fmax=log(huge(y(1))*0.1)                                             3384
19450 do 19451 lam=1,nlam                                                  3384
      f=a0(lam)                                                            3385
19460 do 19461 k=1,nin(lam)                                                3385
      j=ia(k)                                                              3385
      jb=ix(j)                                                             3385
      je=ix(j+1)-1                                                         3386
      f(jx(jb:je))=f(jx(jb:je))+ca(k,lam)*x(jb:je)                         3387
19461 continue                                                             3388
19462 continue                                                             3388
      f=f+g                                                                3389
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                   3390
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                3391
19451 continue                                                             3392
19452 continue                                                             3392
12320 continue                                                             3392
      deallocate(w,f)                                                      3393
      return                                                               3394
      end                                                                  3395
      subroutine multelnet  (parm,no,ni,nr,x,y,w,jd,vp,cl,ne,nx,nlam,flm   3398 
     *in,ulam,thr,isd,jsd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr
     *)
      implicit double precision(a-h,o-z)                                   3399
      double precision x(no,ni),y(no,nr),w(no),vp(ni),ca(nx,nr,nlam)       3400
      double precision ulam(nlam),a0(nr,nlam),rsq(nlam),alm(nlam),cl(2,n   3401 
     *i)
      integer jd(*),ia(nx),nin(nlam)                                       3402
      double precision, dimension (:), allocatable :: vq;                       
      if(maxval(vp) .gt. 0.0)goto 19481                                    3405
      jerr=10000                                                           3405
      return                                                               3405
19481 continue                                                             3406
      allocate(vq(1:ni),stat=jerr)                                         3406
      if(jerr.ne.0) return                                                 3407
      vq=max(0d0,vp)                                                       3407
      vq=vq*ni/sum(vq)                                                     3408
      call multelnetn(parm,no,ni,nr,x,y,w,jd,vq,cl,ne,nx,nlam,flmin,ulam   3410 
     *,thr,isd,  jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      deallocate(vq)                                                       3411
      return                                                               3412
      end                                                                  3413
      subroutine multelnetn (parm,no,ni,nr,x,y,w,jd,vp,cl,ne,nx,nlam,flm   3415 
     *in,ulam,thr,  isd,jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr
     *)
      implicit double precision(a-h,o-z)                                   3416
      double precision vp(ni),x(no,ni),y(no,nr),w(no),ulam(nlam),cl(2,ni   3417 
     *)
      double precision ca(nx,nr,nlam),a0(nr,nlam),rsq(nlam),alm(nlam)      3418
      integer jd(*),ia(nx),nin(nlam)                                       3419
      double precision, dimension (:), allocatable :: xm,xs,xv,ym,ys            
      integer, dimension (:), allocatable :: ju                                 
      double precision, dimension (:,:,:), allocatable :: clt                   
      allocate(clt(1:2,1:nr,1:ni),stat=jerr);                                   
      if(jerr.ne.0) return                                                      
      allocate(xm(1:ni),stat=jerr)                                         3427
      if(jerr.ne.0) return                                                 3428
      allocate(xs(1:ni),stat=jerr)                                         3429
      if(jerr.ne.0) return                                                 3430
      allocate(ym(1:nr),stat=jerr)                                         3431
      if(jerr.ne.0) return                                                 3432
      allocate(ys(1:nr),stat=jerr)                                         3433
      if(jerr.ne.0) return                                                 3434
      allocate(ju(1:ni),stat=jerr)                                         3435
      if(jerr.ne.0) return                                                 3436
      allocate(xv(1:ni),stat=jerr)                                         3437
      if(jerr.ne.0) return                                                 3438
      call chkvars(no,ni,x,ju)                                             3439
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 3440
      if(maxval(ju) .gt. 0)goto 19501                                      3440
      jerr=7777                                                            3440
      return                                                               3440
19501 continue                                                             3441
      call multstandard1(no,ni,nr,x,y,w,isd,jsd,intr,ju,xm,xs,ym,ys,xv,y   3442 
     *s0,jerr)
      if(jerr.ne.0) return                                                 3443
19510 do 19511 j=1,ni                                                      3443
19520 do 19521 k=1,nr                                                      3443
19530 do 19531 i=1,2                                                       3443
      clt(i,k,j)=cl(i,j)                                                   3443
19531 continue                                                             3443
19532 continue                                                             3443
19521 continue                                                             3443
19522 continue                                                             3443
19511 continue                                                             3444
19512 continue                                                             3444
      if(isd .le. 0)goto 19551                                             3444
19560 do 19561 j=1,ni                                                      3444
19570 do 19571 k=1,nr                                                      3444
19580 do 19581 i=1,2                                                       3444
      clt(i,k,j)=clt(i,k,j)*xs(j)                                          3444
19581 continue                                                             3444
19582 continue                                                             3444
19571 continue                                                             3444
19572 continue                                                             3444
19561 continue                                                             3444
19562 continue                                                             3444
19551 continue                                                             3445
      if(jsd .le. 0)goto 19601                                             3445
19610 do 19611 j=1,ni                                                      3445
19620 do 19621 k=1,nr                                                      3445
19630 do 19631 i=1,2                                                       3445
      clt(i,k,j)=clt(i,k,j)/ys(k)                                          3445
19631 continue                                                             3445
19632 continue                                                             3445
19621 continue                                                             3445
19622 continue                                                             3445
19611 continue                                                             3445
19612 continue                                                             3445
19601 continue                                                             3446
      call multelnet2(parm,ni,nr,ju,vp,clt,y,no,ne,nx,x,nlam,flmin,ulam,   3448 
     *thr,maxit,xv,  ys0,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 3449
19640 do 19641 k=1,lmu                                                     3449
      nk=nin(k)                                                            3450
19650 do 19651 j=1,nr                                                      3451
19660 do 19661 l=1,nk                                                      3451
      ca(l,j,k)=ys(j)*ca(l,j,k)/xs(ia(l))                                  3451
19661 continue                                                             3452
19662 continue                                                             3452
      if(intr .ne. 0)goto 19681                                            3452
      a0(j,k)=0.0                                                          3452
      goto 19691                                                           3453
19681 continue                                                             3453
      a0(j,k)=ym(j)-dot_product(ca(1:nk,j,k),xm(ia(1:nk)))                 3453
19691 continue                                                             3454
19671 continue                                                             3454
19651 continue                                                             3455
19652 continue                                                             3455
19641 continue                                                             3456
19642 continue                                                             3456
      deallocate(xm,xs,ym,ys,ju,xv,clt)                                    3457
      return                                                               3458
      end                                                                  3459
      subroutine multstandard1  (no,ni,nr,x,y,w,isd,jsd,intr,ju,xm,xs,ym   3461 
     *,ys,xv,ys0,jerr)
      implicit double precision(a-h,o-z)                                   3462
      double precision x(no,ni),y(no,nr),w(no),xm(ni),xs(ni),xv(ni),ym(n   3463 
     *r),ys(nr)
      integer ju(ni)                                                       3464
      double precision, dimension (:), allocatable :: v                         
      allocate(v(1:no),stat=jerr)                                          3467
      if(jerr.ne.0) return                                                 3468
      w=w/sum(w)                                                           3468
      v=sqrt(w)                                                            3469
      if(intr .ne. 0)goto 19711                                            3470
19720 do 19721 j=1,ni                                                      3470
      if(ju(j).eq.0)goto 19721                                             3470
      xm(j)=0.0                                                            3470
      x(:,j)=v*x(:,j)                                                      3471
      z=dot_product(x(:,j),x(:,j))                                         3472
      if(isd .le. 0)goto 19741                                             3472
      xbq=dot_product(v,x(:,j))**2                                         3472
      vc=z-xbq                                                             3473
      xs(j)=sqrt(vc)                                                       3473
      x(:,j)=x(:,j)/xs(j)                                                  3473
      xv(j)=1.0+xbq/vc                                                     3474
      goto 19751                                                           3475
19741 continue                                                             3475
      xs(j)=1.0                                                            3475
      xv(j)=z                                                              3475
19751 continue                                                             3476
19731 continue                                                             3476
19721 continue                                                             3477
19722 continue                                                             3477
      ys0=0.0                                                              3478
19760 do 19761 j=1,nr                                                      3478
      ym(j)=0.0                                                            3478
      y(:,j)=v*y(:,j)                                                      3479
      z=dot_product(y(:,j),y(:,j))                                         3480
      if(jsd .le. 0)goto 19781                                             3480
      u=z-dot_product(v,y(:,j))**2                                         3480
      ys0=ys0+z/u                                                          3481
      ys(j)=sqrt(u)                                                        3481
      y(:,j)=y(:,j)/ys(j)                                                  3482
      goto 19791                                                           3483
19781 continue                                                             3483
      ys(j)=1.0                                                            3483
      ys0=ys0+z                                                            3483
19791 continue                                                             3484
19771 continue                                                             3484
19761 continue                                                             3485
19762 continue                                                             3485
      go to 10720                                                          3486
19711 continue                                                             3487
19800 do 19801 j=1,ni                                                      3487
      if(ju(j).eq.0)goto 19801                                             3488
      xm(j)=dot_product(w,x(:,j))                                          3488
      x(:,j)=v*(x(:,j)-xm(j))                                              3489
      xv(j)=dot_product(x(:,j),x(:,j))                                     3489
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       3490
19801 continue                                                             3491
19802 continue                                                             3491
      if(isd .ne. 0)goto 19821                                             3491
      xs=1.0                                                               3491
      goto 19831                                                           3492
19821 continue                                                             3492
19840 do 19841 j=1,ni                                                      3492
      if(ju(j).eq.0)goto 19841                                             3492
      x(:,j)=x(:,j)/xs(j)                                                  3492
19841 continue                                                             3493
19842 continue                                                             3493
      xv=1.0                                                               3494
19831 continue                                                             3495
19811 continue                                                             3495
      ys0=0.0                                                              3496
19850 do 19851 j=1,nr                                                      3497
      ym(j)=dot_product(w,y(:,j))                                          3497
      y(:,j)=v*(y(:,j)-ym(j))                                              3498
      z=dot_product(y(:,j),y(:,j))                                         3499
      if(jsd .le. 0)goto 19871                                             3499
      ys(j)=sqrt(z)                                                        3499
      y(:,j)=y(:,j)/ys(j)                                                  3499
      goto 19881                                                           3500
19871 continue                                                             3500
      ys0=ys0+z                                                            3500
19881 continue                                                             3501
19861 continue                                                             3501
19851 continue                                                             3502
19852 continue                                                             3502
      if(jsd .ne. 0)goto 19901                                             3502
      ys=1.0                                                               3502
      goto 19911                                                           3502
19901 continue                                                             3502
      ys0=nr                                                               3502
19911 continue                                                             3503
19891 continue                                                             3503
10720 continue                                                             3503
      deallocate(v)                                                        3504
      return                                                               3505
      end                                                                  3506
      subroutine multelnet2(beta,ni,nr,ju,vp,cl,y,no,ne,nx,x,nlam,flmin,   3508 
     *ulam,thri,  maxit,xv,ys0,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      implicit double precision(a-h,o-z)                                   3509
      double precision vp(ni),y(no,nr),x(no,ni),ulam(nlam),ao(nx,nr,nlam   3510 
     *)
      double precision rsqo(nlam),almo(nlam),xv(ni),cl(2,nr,ni)            3511
      integer ju(ni),ia(nx),kin(nlam)                                      3512
      double precision, dimension (:), allocatable :: g,gk,del,gj               
      integer, dimension (:), allocatable :: mm,ix,isc                          
      double precision, dimension (:,:), allocatable :: a                       
      allocate(a(1:nr,1:ni),stat=jerr)                                          
      if(jerr.ne.0) return                                                      
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)               3520
      allocate(gj(1:nr),stat=jerr)                                         3521
      if(jerr.ne.0) return                                                 3522
      allocate(gk(1:nr),stat=jerr)                                         3523
      if(jerr.ne.0) return                                                 3524
      allocate(del(1:nr),stat=jerr)                                        3525
      if(jerr.ne.0) return                                                 3526
      allocate(mm(1:ni),stat=jerr)                                         3527
      if(jerr.ne.0) return                                                 3528
      allocate(g(1:ni),stat=jerr)                                          3529
      if(jerr.ne.0) return                                                 3530
      allocate(ix(1:ni),stat=jerr)                                         3531
      if(jerr.ne.0) return                                                 3532
      allocate(isc(1:nr),stat=jerr)                                        3533
      if(jerr.ne.0) return                                                 3534
      bta=beta                                                             3534
      omb=1.0-bta                                                          3534
      ix=0                                                                 3534
      thr=thri*ys0/nr                                                      3536
      alf=1.0                                                              3538
      if(flmin .ge. 1.0)goto 19931                                         3538
      eqs=max(eps,flmin)                                                   3538
      alf=eqs**(1.0/(nlam-1))                                              3538
19931 continue                                                             3539
      rsq=ys0                                                              3539
      a=0.0                                                                3539
      mm=0                                                                 3539
      nlp=0                                                                3539
      nin=nlp                                                              3539
      iz=0                                                                 3539
      mnl=min(mnlam,nlam)                                                  3539
      alm=0.0                                                              3540
19940 do 19941 j=1,ni                                                      3540
      if(ju(j).eq.0)goto 19941                                             3540
      g(j)=0.0                                                             3541
19950 do 19951 k=1,nr                                                      3541
      g(j)=g(j)+dot_product(y(:,k),x(:,j))**2                              3541
19951 continue                                                             3542
19952 continue                                                             3542
      g(j)=sqrt(g(j))                                                      3543
19941 continue                                                             3544
19942 continue                                                             3544
19960 do 19961 m=1,nlam                                                    3544
      alm0=alm                                                             3545
      if(flmin .lt. 1.0)goto 19981                                         3545
      alm=ulam(m)                                                          3545
      goto 19971                                                           3546
19981 if(m .le. 2)goto 19991                                               3546
      alm=alm*alf                                                          3546
      goto 19971                                                           3547
19991 if(m .ne. 1)goto 20001                                               3547
      alm=big                                                              3547
      goto 20011                                                           3548
20001 continue                                                             3548
      alm0=0.0                                                             3549
20020 do 20021 j=1,ni                                                      3549
      if(ju(j).eq.0)goto 20021                                             3550
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           3551
20021 continue                                                             3552
20022 continue                                                             3552
      alm0=alm0/max(bta,1.0d-3)                                            3552
      alm=alf*alm0                                                         3553
20011 continue                                                             3554
19971 continue                                                             3554
      dem=alm*omb                                                          3554
      ab=alm*bta                                                           3554
      rsq0=rsq                                                             3554
      jz=1                                                                 3555
      tlam=bta*(2.0*alm-alm0)                                              3556
20030 do 20031 k=1,ni                                                      3556
      if(ix(k).eq.1)goto 20031                                             3556
      if(ju(k).eq.0)goto 20031                                             3557
      if(g(k).gt.tlam*vp(k)) ix(k)=1                                       3558
20031 continue                                                             3559
20032 continue                                                             3559
20040 continue                                                             3559
20041 continue                                                             3559
      if(iz*jz.ne.0) go to 10360                                           3560
11020 continue                                                             3560
      nlp=nlp+1                                                            3560
      dlx=0.0                                                              3561
20050 do 20051 k=1,ni                                                      3561
      if(ix(k).eq.0)goto 20051                                             3561
      gkn=0.0                                                              3562
20060 do 20061 j=1,nr                                                      3562
      gj(j)=dot_product(y(:,j),x(:,k))                                     3563
      gk(j)=gj(j)+a(j,k)*xv(k)                                             3563
      gkn=gkn+gk(j)**2                                                     3565
20061 continue                                                             3565
20062 continue                                                             3565
      gkn=sqrt(gkn)                                                        3565
      u=1.0-ab*vp(k)/gkn                                                   3565
      del=a(:,k)                                                           3566
      if(u .gt. 0.0)goto 20081                                             3566
      a(:,k)=0.0                                                           3566
      goto 20091                                                           3567
20081 continue                                                             3567
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      3568
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)   3570 
     *,isc,jerr)
      if(jerr.ne.0) return                                                 3571
20091 continue                                                             3572
20071 continue                                                             3572
      del=a(:,k)-del                                                       3572
      if(maxval(abs(del)).le.0.0)goto 20051                                3573
20100 do 20101 j=1,nr                                                      3573
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              3574
      y(:,j)=y(:,j)-del(j)*x(:,k)                                          3574
      dlx=max(dlx,xv(k)*del(j)**2)                                         3575
20101 continue                                                             3576
20102 continue                                                             3576
      if(mm(k) .ne. 0)goto 20121                                           3576
      nin=nin+1                                                            3576
      if(nin.gt.nx)goto 20052                                              3577
      mm(k)=nin                                                            3577
      ia(nin)=k                                                            3578
20121 continue                                                             3579
20051 continue                                                             3580
20052 continue                                                             3580
      if(nin.gt.nx)goto 20042                                              3581
      if(dlx .ge. thr)goto 20141                                           3581
      ixx=0                                                                3582
20150 do 20151 k=1,ni                                                      3582
      if(ix(k).eq.1)goto 20151                                             3582
      if(ju(k).eq.0)goto 20151                                             3582
      g(k)=0.0                                                             3583
20160 do 20161 j=1,nr                                                      3583
      g(k)=g(k)+dot_product(y(:,j),x(:,k))**2                              3583
20161 continue                                                             3584
20162 continue                                                             3584
      g(k)=sqrt(g(k))                                                      3585
      if(g(k) .le. ab*vp(k))goto 20181                                     3585
      ix(k)=1                                                              3585
      ixx=1                                                                3585
20181 continue                                                             3586
20151 continue                                                             3587
20152 continue                                                             3587
      if(ixx.eq.1) go to 11020                                             3588
      goto 20042                                                           3589
20141 continue                                                             3590
      if(nlp .le. maxit)goto 20201                                         3590
      jerr=-m                                                              3590
      return                                                               3590
20201 continue                                                             3591
10360 continue                                                             3591
      iz=1                                                                 3592
20210 continue                                                             3592
20211 continue                                                             3592
      nlp=nlp+1                                                            3592
      dlx=0.0                                                              3593
20220 do 20221 l=1,nin                                                     3593
      k=ia(l)                                                              3593
      gkn=0.0                                                              3594
20230 do 20231 j=1,nr                                                      3594
      gj(j)=dot_product(y(:,j),x(:,k))                                     3595
      gk(j)=gj(j)+a(j,k)*xv(k)                                             3595
      gkn=gkn+gk(j)**2                                                     3597
20231 continue                                                             3597
20232 continue                                                             3597
      gkn=sqrt(gkn)                                                        3597
      u=1.0-ab*vp(k)/gkn                                                   3597
      del=a(:,k)                                                           3598
      if(u .gt. 0.0)goto 20251                                             3598
      a(:,k)=0.0                                                           3598
      goto 20261                                                           3599
20251 continue                                                             3599
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      3600
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)   3602 
     *,isc,jerr)
      if(jerr.ne.0) return                                                 3603
20261 continue                                                             3604
20241 continue                                                             3604
      del=a(:,k)-del                                                       3604
      if(maxval(abs(del)).le.0.0)goto 20221                                3605
20270 do 20271 j=1,nr                                                      3605
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              3606
      y(:,j)=y(:,j)-del(j)*x(:,k)                                          3606
      dlx=max(dlx,xv(k)*del(j)**2)                                         3607
20271 continue                                                             3608
20272 continue                                                             3608
20221 continue                                                             3609
20222 continue                                                             3609
      if(dlx.lt.thr)goto 20212                                             3609
      if(nlp .le. maxit)goto 20291                                         3609
      jerr=-m                                                              3609
      return                                                               3609
20291 continue                                                             3610
      goto 20211                                                           3611
20212 continue                                                             3611
      jz=0                                                                 3612
      goto 20041                                                           3613
20042 continue                                                             3613
      if(nin .le. nx)goto 20311                                            3613
      jerr=-10000-m                                                        3613
      goto 19962                                                           3613
20311 continue                                                             3614
      if(nin .le. 0)goto 20331                                             3614
20340 do 20341 j=1,nr                                                      3614
      ao(1:nin,j,m)=a(j,ia(1:nin))                                         3614
20341 continue                                                             3614
20342 continue                                                             3614
20331 continue                                                             3615
      kin(m)=nin                                                           3616
      rsqo(m)=1.0-rsq/ys0                                                  3616
      almo(m)=alm                                                          3616
      lmu=m                                                                3617
      if(m.lt.mnl)goto 19961                                               3617
      if(flmin.ge.1.0)goto 19961                                           3618
      me=0                                                                 3618
20350 do 20351 j=1,nin                                                     3618
      if(ao(j,1,m).ne.0.0) me=me+1                                         3618
20351 continue                                                             3618
20352 continue                                                             3618
      if(me.gt.ne)goto 19962                                               3619
      if(rsq0-rsq.lt.sml*rsq)goto 19962                                    3619
      if(rsqo(m).gt.rsqmax)goto 19962                                      3620
19961 continue                                                             3621
19962 continue                                                             3621
      deallocate(a,mm,g,ix,del,gj,gk)                                      3622
      return                                                               3623
      end                                                                  3624
      subroutine chkbnds(nr,gk,gkn,xv,cl,al1,al2,a,isc,jerr)               3625
      implicit double precision(a-h,o-z)                                   3626
      double precision gk(nr),cl(2,nr),a(nr)                               3626
      integer isc(nr)                                                      3627
      kerr=0                                                               3627
      al1p=1.0+al1/xv                                                      3627
      al2p=al2/xv                                                          3627
      isc=0                                                                3628
      gsq=gkn**2                                                           3628
      asq=dot_product(a,a)                                                 3628
      usq=0.0                                                              3630
      u=0.0                                                                3630
      kn=-1                                                                3632
20360 continue                                                             3632
20361 continue                                                             3632
      vmx=0.0                                                              3633
20370 do 20371 k=1,nr                                                      3633
      v=max(a(k)-cl(2,k),cl(1,k)-a(k))                                     3634
      if(v .le. vmx)goto 20391                                             3634
      vmx=v                                                                3634
      kn=k                                                                 3634
20391 continue                                                             3635
20371 continue                                                             3636
20372 continue                                                             3636
      if(vmx.le.0.0)goto 20362                                             3636
      if(isc(kn).ne.0)goto 20362                                           3637
      gsq=gsq-gk(kn)**2                                                    3637
      g=sqrt(gsq)/xv                                                       3638
      if(a(kn).lt.cl(1,kn)) u=cl(1,kn)                                     3638
      if(a(kn).gt.cl(2,kn)) u=cl(2,kn)                                     3639
      usq=usq+u**2                                                         3640
      if(usq .ne. 0.0)goto 20411                                           3640
      b=max(0d0,(g-al2p)/al1p)                                             3640
      goto 20421                                                           3641
20411 continue                                                             3641
      b0=sqrt(asq-a(kn)**2)                                                3642
      b=bnorm(b0,al1p,al2p,g,usq,kerr)                                     3642
      if(kerr.ne.0)goto 20362                                              3643
20421 continue                                                             3644
20401 continue                                                             3644
      asq=usq+b**2                                                         3644
      if(asq .gt. 0.0)goto 20441                                           3644
      a=0.0                                                                3644
      goto 20362                                                           3644
20441 continue                                                             3645
      a(kn)=u                                                              3645
      isc(kn)=1                                                            3645
      f=1.0/(xv*(al1p+al2p/sqrt(asq)))                                     3646
20450 do 20451 j=1,nr                                                      3646
      if(isc(j).eq.0) a(j)=f*gk(j)                                         3646
20451 continue                                                             3647
20452 continue                                                             3647
      goto 20361                                                           3648
20362 continue                                                             3648
      if(kerr.ne.0) jerr=kerr                                              3649
      return                                                               3650
      end                                                                  3651
      subroutine chkbnds1(nr,gk,gkn,xv,cl1,cl2,al1,al2,a,isc,jerr)         3652
      implicit double precision(a-h,o-z)                                   3653
      double precision gk(nr),a(nr)                                        3653
      integer isc(nr)                                                      3654
      kerr=0                                                               3654
      al1p=1.0+al1/xv                                                      3654
      al2p=al2/xv                                                          3654
      isc=0                                                                3655
      gsq=gkn**2                                                           3655
      asq=dot_product(a,a)                                                 3655
      usq=0.0                                                              3657
      u=0.0                                                                3657
      kn=-1                                                                3659
20460 continue                                                             3659
20461 continue                                                             3659
      vmx=0.0                                                              3660
20470 do 20471 k=1,nr                                                      3660
      v=max(a(k)-cl2,cl1-a(k))                                             3661
      if(v .le. vmx)goto 20491                                             3661
      vmx=v                                                                3661
      kn=k                                                                 3661
20491 continue                                                             3662
20471 continue                                                             3663
20472 continue                                                             3663
      if(vmx.le.0.0)goto 20462                                             3663
      if(isc(kn).ne.0)goto 20462                                           3664
      gsq=gsq-gk(kn)**2                                                    3664
      g=sqrt(gsq)/xv                                                       3665
      if(a(kn).lt.cl1) u=cl1                                               3665
      if(a(kn).gt.cl2) u=cl2                                               3666
      usq=usq+u**2                                                         3667
      if(usq .ne. 0.0)goto 20511                                           3667
      b=max(0d0,(g-al2p)/al1p)                                             3667
      goto 20521                                                           3668
20511 continue                                                             3668
      b0=sqrt(asq-a(kn)**2)                                                3669
      b=bnorm(b0,al1p,al2p,g,usq,kerr)                                     3669
      if(kerr.ne.0)goto 20462                                              3670
20521 continue                                                             3671
20501 continue                                                             3671
      asq=usq+b**2                                                         3671
      if(asq .gt. 0.0)goto 20541                                           3671
      a=0.0                                                                3671
      goto 20462                                                           3671
20541 continue                                                             3672
      a(kn)=u                                                              3672
      isc(kn)=1                                                            3672
      f=1.0/(xv*(al1p+al2p/sqrt(asq)))                                     3673
20550 do 20551 j=1,nr                                                      3673
      if(isc(j).eq.0) a(j)=f*gk(j)                                         3673
20551 continue                                                             3674
20552 continue                                                             3674
      goto 20461                                                           3675
20462 continue                                                             3675
      if(kerr.ne.0) jerr=kerr                                              3676
      return                                                               3677
      end                                                                  3678
      function bnorm(b0,al1p,al2p,g,usq,jerr)                              3679
      implicit double precision(a-h,o-z)                                   3680
      data thr,mxit /1.0d-10,100/                                          3681
      b=b0                                                                 3681
      zsq=b**2+usq                                                         3681
      if(zsq .gt. 0.0)goto 20571                                           3681
      bnorm=0.0                                                            3681
      return                                                               3681
20571 continue                                                             3682
      z=sqrt(zsq)                                                          3682
      f=b*(al1p+al2p/z)-g                                                  3682
      jerr=0                                                               3683
20580 do 20581 it=1,mxit                                                   3683
      b=b-f/(al1p+al2p*usq/(z*zsq))                                        3684
      zsq=b**2+usq                                                         3684
      if(zsq .gt. 0.0)goto 20601                                           3684
      bnorm=0.0                                                            3684
      return                                                               3684
20601 continue                                                             3685
      z=sqrt(zsq)                                                          3685
      f=b*(al1p+al2p/z)-g                                                  3686
      if(abs(f).le.thr)goto 20582                                          3686
      if(b .gt. 0.0)goto 20621                                             3686
      b=0.0                                                                3686
      goto 20582                                                           3686
20621 continue                                                             3687
20581 continue                                                             3688
20582 continue                                                             3688
      bnorm=b                                                              3688
      if(it.ge.mxit) jerr=90000                                            3689
      return                                                               3691
      entry chg_bnorm(arg,irg)                                             3691
      bnorm = 0.0                                                          3691
      thr=arg                                                              3691
      mxit=irg                                                             3691
      return                                                               3692
      entry get_bnorm(arg,irg)                                             3692
      bnorm = 0.0                                                          3692
      arg=thr                                                              3692
      irg=mxit                                                             3692
      return                                                               3694
      end                                                                  3695
      subroutine multsolns(ni,nx,nr,lmu,a,ia,nin,b)                        3696
      implicit double precision(a-h,o-z)                                   3697
      double precision a(nx,nr,lmu),b(ni,nr,lmu)                           3697
      integer ia(nx),nin(lmu)                                              3698
20630 do 20631 lam=1,lmu                                                   3698
      call multuncomp(ni,nr,nx,a(1,1,lam),ia,nin(lam),b(1,1,lam))          3698
20631 continue                                                             3699
20632 continue                                                             3699
      return                                                               3700
      end                                                                  3701
      subroutine multuncomp(ni,nr,nx,ca,ia,nin,a)                          3702
      implicit double precision(a-h,o-z)                                   3703
      double precision ca(nx,nr),a(ni,nr)                                  3703
      integer ia(nx)                                                       3704
      a=0.0                                                                3705
      if(nin .le. 0)goto 20651                                             3705
20660 do 20661 j=1,nr                                                      3705
      a(ia(1:nin),j)=ca(1:nin,j)                                           3705
20661 continue                                                             3705
20662 continue                                                             3705
20651 continue                                                             3706
      return                                                               3707
      end                                                                  3708
      subroutine multmodval(nx,nr,a0,ca,ia,nin,n,x,f)                      3709
      implicit double precision(a-h,o-z)                                   3710
      double precision a0(nr),ca(nx,nr),x(n,*),f(nr,n)                     3710
      integer ia(nx)                                                       3711
20670 do 20671 i=1,n                                                       3711
      f(:,i)=a0                                                            3711
20671 continue                                                             3711
20672 continue                                                             3711
      if(nin.le.0) return                                                  3712
20680 do 20681 i=1,n                                                       3712
20690 do 20691 j=1,nr                                                      3712
      f(j,i)=f(j,i)+dot_product(ca(1:nin,j),x(i,ia(1:nin)))                3712
20691 continue                                                             3712
20692 continue                                                             3712
20681 continue                                                             3713
20682 continue                                                             3713
      return                                                               3714
      end                                                                  3715
      subroutine multspelnet  (parm,no,ni,nr,x,ix,jx,y,w,jd,vp,cl,ne,nx,   3718 
     *nlam,flmin,ulam,thr,isd,  jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,
     *nlp,jerr)
      implicit double precision(a-h,o-z)                                   3719
      double precision x(*),y(no,nr),w(no),vp(ni),ulam(nlam),cl(2,ni)      3720
      double precision ca(nx,nr,nlam),a0(nr,nlam),rsq(nlam),alm(nlam)      3721
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           3722
      double precision, dimension (:), allocatable :: vq;                       
      if(maxval(vp) .gt. 0.0)goto 20711                                    3725
      jerr=10000                                                           3725
      return                                                               3725
20711 continue                                                             3726
      allocate(vq(1:ni),stat=jerr)                                         3726
      if(jerr.ne.0) return                                                 3727
      vq=max(0d0,vp)                                                       3727
      vq=vq*ni/sum(vq)                                                     3728
      call multspelnetn(parm,no,ni,nr,x,ix,jx,y,w,jd,vq,cl,ne,nx,nlam,fl   3730 
     *min,  ulam,thr,isd,jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jer
     *r)
      deallocate(vq)                                                       3731
      return                                                               3732
      end                                                                  3733
      subroutine multspelnetn(parm,no,ni,nr,x,ix,jx,y,w,jd,vp,cl,ne,nx,n   3735 
     *lam,flmin,  ulam,thr,isd,jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,n
     *lp,jerr)
      implicit double precision(a-h,o-z)                                   3736
      double precision x(*),vp(ni),y(no,nr),w(no),ulam(nlam),cl(2,ni)      3737
      double precision ca(nx,nr,nlam),a0(nr,nlam),rsq(nlam),alm(nlam)      3738
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           3739
      double precision, dimension (:), allocatable :: xm,xs,xv,ym,ys            
      integer, dimension (:), allocatable :: ju                                 
      double precision, dimension (:,:,:), allocatable :: clt                   
      allocate(clt(1:2,1:nr,1:ni),stat=jerr)                                    
      if(jerr.ne.0) return                                                      
      allocate(xm(1:ni),stat=jerr)                                         3747
      if(jerr.ne.0) return                                                 3748
      allocate(xs(1:ni),stat=jerr)                                         3749
      if(jerr.ne.0) return                                                 3750
      allocate(ym(1:nr),stat=jerr)                                         3751
      if(jerr.ne.0) return                                                 3752
      allocate(ys(1:nr),stat=jerr)                                         3753
      if(jerr.ne.0) return                                                 3754
      allocate(ju(1:ni),stat=jerr)                                         3755
      if(jerr.ne.0) return                                                 3756
      allocate(xv(1:ni),stat=jerr)                                         3757
      if(jerr.ne.0) return                                                 3758
      call spchkvars(no,ni,x,ix,ju)                                        3759
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 3760
      if(maxval(ju) .gt. 0)goto 20731                                      3760
      jerr=7777                                                            3760
      return                                                               3760
20731 continue                                                             3761
      call multspstandard1(no,ni,nr,x,ix,jx,y,w,ju,isd,jsd,intr,  xm,xs,   3763 
     *ym,ys,xv,ys0,jerr)
      if(jerr.ne.0) return                                                 3764
20740 do 20741 j=1,ni                                                      3764
20750 do 20751 k=1,nr                                                      3764
20760 do 20761 i=1,2                                                       3764
      clt(i,k,j)=cl(i,j)                                                   3764
20761 continue                                                             3764
20762 continue                                                             3764
20751 continue                                                             3764
20752 continue                                                             3764
20741 continue                                                             3765
20742 continue                                                             3765
      if(isd .le. 0)goto 20781                                             3765
20790 do 20791 j=1,ni                                                      3765
20800 do 20801 k=1,nr                                                      3765
20810 do 20811 i=1,2                                                       3765
      clt(i,k,j)=clt(i,k,j)*xs(j)                                          3765
20811 continue                                                             3765
20812 continue                                                             3765
20801 continue                                                             3765
20802 continue                                                             3765
20791 continue                                                             3765
20792 continue                                                             3765
20781 continue                                                             3766
      if(jsd .le. 0)goto 20831                                             3766
20840 do 20841 j=1,ni                                                      3766
20850 do 20851 k=1,nr                                                      3766
20860 do 20861 i=1,2                                                       3766
      clt(i,k,j)=clt(i,k,j)/ys(k)                                          3766
20861 continue                                                             3766
20862 continue                                                             3766
20851 continue                                                             3766
20852 continue                                                             3766
20841 continue                                                             3766
20842 continue                                                             3766
20831 continue                                                             3767
      call multspelnet2(parm,ni,nr,y,w,no,ne,nx,x,ix,jx,ju,vp,clt,nlam,f   3769 
     *lmin,  ulam,thr,maxit,xm,xs,xv,ys0,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 3770
20870 do 20871 k=1,lmu                                                     3770
      nk=nin(k)                                                            3771
20880 do 20881 j=1,nr                                                      3772
20890 do 20891 l=1,nk                                                      3772
      ca(l,j,k)=ys(j)*ca(l,j,k)/xs(ia(l))                                  3772
20891 continue                                                             3773
20892 continue                                                             3773
      if(intr .ne. 0)goto 20911                                            3773
      a0(j,k)=0.0                                                          3773
      goto 20921                                                           3774
20911 continue                                                             3774
      a0(j,k)=ym(j)-dot_product(ca(1:nk,j,k),xm(ia(1:nk)))                 3774
20921 continue                                                             3775
20901 continue                                                             3775
20881 continue                                                             3776
20882 continue                                                             3776
20871 continue                                                             3777
20872 continue                                                             3777
      deallocate(xm,xs,ym,ys,ju,xv,clt)                                    3778
      return                                                               3779
      end                                                                  3780
      subroutine multspstandard1(no,ni,nr,x,ix,jx,y,w,ju,isd,jsd,intr,     3782 
     *xm,xs,ym,ys,xv,ys0,jerr)
      implicit double precision(a-h,o-z)                                   3783
      double precision x(*),y(no,nr),w(no),xm(ni),xs(ni),xv(ni),ym(nr),y   3784 
     *s(nr)
      integer ix(*),jx(*),ju(ni)                                           3785
      w=w/sum(w)                                                           3786
      if(intr .ne. 0)goto 20941                                            3787
20950 do 20951 j=1,ni                                                      3787
      if(ju(j).eq.0)goto 20951                                             3787
      xm(j)=0.0                                                            3787
      jb=ix(j)                                                             3787
      je=ix(j+1)-1                                                         3788
      z=dot_product(w(jx(jb:je)),x(jb:je)**2)                              3789
      if(isd .le. 0)goto 20971                                             3789
      xbq=dot_product(w(jx(jb:je)),x(jb:je))**2                            3789
      vc=z-xbq                                                             3790
      xs(j)=sqrt(vc)                                                       3790
      xv(j)=1.0+xbq/vc                                                     3791
      goto 20981                                                           3792
20971 continue                                                             3792
      xs(j)=1.0                                                            3792
      xv(j)=z                                                              3792
20981 continue                                                             3793
20961 continue                                                             3793
20951 continue                                                             3794
20952 continue                                                             3794
      ys0=0.0                                                              3795
20990 do 20991 j=1,nr                                                      3795
      ym(j)=0.0                                                            3795
      z=dot_product(w,y(:,j)**2)                                           3796
      if(jsd .le. 0)goto 21011                                             3796
      u=z-dot_product(w,y(:,j))**2                                         3796
      ys0=ys0+z/u                                                          3797
      ys(j)=sqrt(u)                                                        3797
      y(:,j)=y(:,j)/ys(j)                                                  3798
      goto 21021                                                           3799
21011 continue                                                             3799
      ys(j)=1.0                                                            3799
      ys0=ys0+z                                                            3799
21021 continue                                                             3800
21001 continue                                                             3800
20991 continue                                                             3801
20992 continue                                                             3801
      return                                                               3802
20941 continue                                                             3803
21030 do 21031 j=1,ni                                                      3803
      if(ju(j).eq.0)goto 21031                                             3804
      jb=ix(j)                                                             3804
      je=ix(j+1)-1                                                         3804
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             3805
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 3806
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       3807
21031 continue                                                             3808
21032 continue                                                             3808
      if(isd .ne. 0)goto 21051                                             3808
      xs=1.0                                                               3808
      goto 21061                                                           3808
21051 continue                                                             3808
      xv=1.0                                                               3808
21061 continue                                                             3809
21041 continue                                                             3809
      ys0=0.0                                                              3810
21070 do 21071 j=1,nr                                                      3811
      ym(j)=dot_product(w,y(:,j))                                          3811
      y(:,j)=y(:,j)-ym(j)                                                  3812
      z=dot_product(w,y(:,j)**2)                                           3813
      if(jsd .le. 0)goto 21091                                             3813
      ys(j)=sqrt(z)                                                        3813
      y(:,j)=y(:,j)/ys(j)                                                  3813
      goto 21101                                                           3814
21091 continue                                                             3814
      ys0=ys0+z                                                            3814
21101 continue                                                             3815
21081 continue                                                             3815
21071 continue                                                             3816
21072 continue                                                             3816
      if(jsd .ne. 0)goto 21121                                             3816
      ys=1.0                                                               3816
      goto 21131                                                           3816
21121 continue                                                             3816
      ys0=nr                                                               3816
21131 continue                                                             3817
21111 continue                                                             3817
      return                                                               3818
      end                                                                  3819
      subroutine multspelnet2(beta,ni,nr,y,w,no,ne,nx,x,ix,jx,ju,vp,cl,n   3821 
     *lam,flmin,  ulam,thri,maxit,xm,xs,xv,ys0,lmu,ao,ia,kin,rsqo,almo,n
     *lp,jerr)
      implicit double precision(a-h,o-z)                                   3822
      double precision y(no,nr),w(no),x(*),vp(ni),ulam(nlam),cl(2,nr,ni)   3823
      double precision ao(nx,nr,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni   3824 
     *),xv(ni)
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                          3825
      double precision, dimension (:), allocatable :: g,gj,gk,del,o             
      integer, dimension (:), allocatable :: mm,iy,isc                          
      double precision, dimension (:,:), allocatable :: a                       
      allocate(a(1:nr,1:ni),stat=jerr)                                          
      if(jerr.ne.0) return                                                      
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)               3833
      allocate(mm(1:ni),stat=jerr)                                         3834
      if(jerr.ne.0) return                                                 3835
      allocate(g(1:ni),stat=jerr)                                          3836
      if(jerr.ne.0) return                                                 3837
      allocate(gj(1:nr),stat=jerr)                                         3838
      if(jerr.ne.0) return                                                 3839
      allocate(gk(1:nr),stat=jerr)                                         3840
      if(jerr.ne.0) return                                                 3841
      allocate(del(1:nr),stat=jerr)                                        3842
      if(jerr.ne.0) return                                                 3843
      allocate(o(1:nr),stat=jerr)                                          3844
      if(jerr.ne.0) return                                                 3845
      allocate(iy(1:ni),stat=jerr)                                         3846
      if(jerr.ne.0) return                                                 3847
      allocate(isc(1:nr),stat=jerr)                                        3848
      if(jerr.ne.0) return                                                 3849
      bta=beta                                                             3849
      omb=1.0-bta                                                          3849
      alm=0.0                                                              3849
      iy=0                                                                 3849
      thr=thri*ys0/nr                                                      3851
      alf=1.0                                                              3853
      if(flmin .ge. 1.0)goto 21151                                         3853
      eqs=max(eps,flmin)                                                   3853
      alf=eqs**(1.0/(nlam-1))                                              3853
21151 continue                                                             3854
      rsq=ys0                                                              3854
      a=0.0                                                                3854
      mm=0                                                                 3854
      o=0.0                                                                3854
      nlp=0                                                                3854
      nin=nlp                                                              3854
      iz=0                                                                 3854
      mnl=min(mnlam,nlam)                                                  3855
21160 do 21161 j=1,ni                                                      3855
      if(ju(j).eq.0)goto 21161                                             3855
      jb=ix(j)                                                             3855
      je=ix(j+1)-1                                                         3855
      g(j)=0.0                                                             3856
21170 do 21171 k=1,nr                                                      3857
      g(j)=g(j)+(dot_product(y(jx(jb:je),k),w(jx(jb:je))*x(jb:je))/xs(j)   3858 
     *)**2
21171 continue                                                             3859
21172 continue                                                             3859
      g(j)=sqrt(g(j))                                                      3860
21161 continue                                                             3861
21162 continue                                                             3861
21180 do 21181 m=1,nlam                                                    3861
      alm0=alm                                                             3862
      if(flmin .lt. 1.0)goto 21201                                         3862
      alm=ulam(m)                                                          3862
      goto 21191                                                           3863
21201 if(m .le. 2)goto 21211                                               3863
      alm=alm*alf                                                          3863
      goto 21191                                                           3864
21211 if(m .ne. 1)goto 21221                                               3864
      alm=big                                                              3864
      goto 21231                                                           3865
21221 continue                                                             3865
      alm0=0.0                                                             3866
21240 do 21241 j=1,ni                                                      3866
      if(ju(j).eq.0)goto 21241                                             3867
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           3868
21241 continue                                                             3869
21242 continue                                                             3869
      alm0=alm0/max(bta,1.0d-3)                                            3869
      alm=alf*alm0                                                         3870
21231 continue                                                             3871
21191 continue                                                             3871
      dem=alm*omb                                                          3871
      ab=alm*bta                                                           3871
      rsq0=rsq                                                             3871
      jz=1                                                                 3872
      tlam=bta*(2.0*alm-alm0)                                              3873
21250 do 21251 k=1,ni                                                      3873
      if(iy(k).eq.1)goto 21251                                             3873
      if(ju(k).eq.0)goto 21251                                             3874
      if(g(k).gt.tlam*vp(k)) iy(k)=1                                       3875
21251 continue                                                             3876
21252 continue                                                             3876
21260 continue                                                             3876
21261 continue                                                             3876
      if(iz*jz.ne.0) go to 10360                                           3877
11020 continue                                                             3877
      nlp=nlp+1                                                            3877
      dlx=0.0                                                              3878
21270 do 21271 k=1,ni                                                      3878
      if(iy(k).eq.0)goto 21271                                             3878
      jb=ix(k)                                                             3878
      je=ix(k+1)-1                                                         3878
      gkn=0.0                                                              3879
21280 do 21281 j=1,nr                                                      3880
      gj(j)=dot_product(y(jx(jb:je),j)+o(j),w(jx(jb:je))*x(jb:je))/xs(k)   3881
      gk(j)=gj(j)+a(j,k)*xv(k)                                             3881
      gkn=gkn+gk(j)**2                                                     3882
21281 continue                                                             3883
21282 continue                                                             3883
      gkn=sqrt(gkn)                                                        3883
      u=1.0-ab*vp(k)/gkn                                                   3883
      del=a(:,k)                                                           3884
      if(u .gt. 0.0)goto 21301                                             3884
      a(:,k)=0.0                                                           3884
      goto 21311                                                           3885
21301 continue                                                             3885
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      3886
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)   3888 
     *,isc,jerr)
      if(jerr.ne.0) return                                                 3889
21311 continue                                                             3890
21291 continue                                                             3890
      del=a(:,k)-del                                                       3890
      if(maxval(abs(del)).le.0.0)goto 21271                                3891
      if(mm(k) .ne. 0)goto 21331                                           3891
      nin=nin+1                                                            3891
      if(nin.gt.nx)goto 21272                                              3892
      mm(k)=nin                                                            3892
      ia(nin)=k                                                            3893
21331 continue                                                             3894
21340 do 21341 j=1,nr                                                      3894
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              3895
      y(jx(jb:je),j)=y(jx(jb:je),j)-del(j)*x(jb:je)/xs(k)                  3896
      o(j)=o(j)+del(j)*xm(k)/xs(k)                                         3896
      dlx=max(xv(k)*del(j)**2,dlx)                                         3897
21341 continue                                                             3898
21342 continue                                                             3898
21271 continue                                                             3899
21272 continue                                                             3899
      if(nin.gt.nx)goto 21262                                              3900
      if(dlx .ge. thr)goto 21361                                           3900
      ixx=0                                                                3901
21370 do 21371 j=1,ni                                                      3901
      if(iy(j).eq.1)goto 21371                                             3901
      if(ju(j).eq.0)goto 21371                                             3902
      jb=ix(j)                                                             3902
      je=ix(j+1)-1                                                         3902
      g(j)=0.0                                                             3903
21380 do 21381 k=1,nr                                                      3903
      g(j)=g(j)+  (dot_product(y(jx(jb:je),k)+o(k),w(jx(jb:je))*x(jb:je)   3905 
     *)/xs(j))**2
21381 continue                                                             3906
21382 continue                                                             3906
      g(j)=sqrt(g(j))                                                      3907
      if(g(j) .le. ab*vp(j))goto 21401                                     3907
      iy(j)=1                                                              3907
      ixx=1                                                                3907
21401 continue                                                             3908
21371 continue                                                             3909
21372 continue                                                             3909
      if(ixx.eq.1) go to 11020                                             3910
      goto 21262                                                           3911
21361 continue                                                             3912
      if(nlp .le. maxit)goto 21421                                         3912
      jerr=-m                                                              3912
      return                                                               3912
21421 continue                                                             3913
10360 continue                                                             3913
      iz=1                                                                 3914
21430 continue                                                             3914
21431 continue                                                             3914
      nlp=nlp+1                                                            3914
      dlx=0.0                                                              3915
21440 do 21441 l=1,nin                                                     3915
      k=ia(l)                                                              3915
      jb=ix(k)                                                             3915
      je=ix(k+1)-1                                                         3915
      gkn=0.0                                                              3916
21450 do 21451 j=1,nr                                                      3916
      gj(j)=  dot_product(y(jx(jb:je),j)+o(j),w(jx(jb:je))*x(jb:je))/xs(   3918 
     *k)
      gk(j)=gj(j)+a(j,k)*xv(k)                                             3918
      gkn=gkn+gk(j)**2                                                     3919
21451 continue                                                             3920
21452 continue                                                             3920
      gkn=sqrt(gkn)                                                        3920
      u=1.0-ab*vp(k)/gkn                                                   3920
      del=a(:,k)                                                           3921
      if(u .gt. 0.0)goto 21471                                             3921
      a(:,k)=0.0                                                           3921
      goto 21481                                                           3922
21471 continue                                                             3922
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      3923
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)   3925 
     *,isc,jerr)
      if(jerr.ne.0) return                                                 3926
21481 continue                                                             3927
21461 continue                                                             3927
      del=a(:,k)-del                                                       3927
      if(maxval(abs(del)).le.0.0)goto 21441                                3928
21490 do 21491 j=1,nr                                                      3928
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              3929
      y(jx(jb:je),j)=y(jx(jb:je),j)-del(j)*x(jb:je)/xs(k)                  3930
      o(j)=o(j)+del(j)*xm(k)/xs(k)                                         3930
      dlx=max(xv(k)*del(j)**2,dlx)                                         3931
21491 continue                                                             3932
21492 continue                                                             3932
21441 continue                                                             3933
21442 continue                                                             3933
      if(dlx.lt.thr)goto 21432                                             3933
      if(nlp .le. maxit)goto 21511                                         3933
      jerr=-m                                                              3933
      return                                                               3933
21511 continue                                                             3934
      goto 21431                                                           3935
21432 continue                                                             3935
      jz=0                                                                 3936
      goto 21261                                                           3937
21262 continue                                                             3937
      if(nin .le. nx)goto 21531                                            3937
      jerr=-10000-m                                                        3937
      goto 21182                                                           3937
21531 continue                                                             3938
      if(nin .le. 0)goto 21551                                             3938
21560 do 21561 j=1,nr                                                      3938
      ao(1:nin,j,m)=a(j,ia(1:nin))                                         3938
21561 continue                                                             3938
21562 continue                                                             3938
21551 continue                                                             3939
      kin(m)=nin                                                           3940
      rsqo(m)=1.0-rsq/ys0                                                  3940
      almo(m)=alm                                                          3940
      lmu=m                                                                3941
      if(m.lt.mnl)goto 21181                                               3941
      if(flmin.ge.1.0)goto 21181                                           3942
      me=0                                                                 3942
21570 do 21571 j=1,nin                                                     3942
      if(ao(j,1,m).ne.0.0) me=me+1                                         3942
21571 continue                                                             3942
21572 continue                                                             3942
      if(me.gt.ne)goto 21182                                               3943
      if(rsq0-rsq.lt.sml*rsq)goto 21182                                    3943
      if(rsqo(m).gt.rsqmax)goto 21182                                      3944
21181 continue                                                             3945
21182 continue                                                             3945
      deallocate(a,mm,g,iy,gj,gk,del,o)                                    3946
      return                                                               3947
      end                                                                  3948
      subroutine multlognetn(parm,no,ni,nc,x,y,g,w,ju,vp,cl,ne,nx,nlam,f   3950 
     *lmin,ulam,  shri,intr,maxit,xv,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jer
     *r)
      implicit double precision(a-h,o-z)                                   3951
      double precision x(no,ni),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam   3952 
     *),cl(2,ni)
      double precision a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),xv(   3953 
     *ni)
      integer ju(ni),m(nx),kin(nlam)                                       3954
      double precision, dimension (:,:), allocatable :: q,r,b,bs                
      double precision, dimension (:), allocatable :: sxp,sxpl,ga,gk,del        
      integer, dimension (:), allocatable :: mm,is,ixx,isc                      
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      if(jerr.ne.0) return                                                      
      allocate(bs(0:ni,1:nc),stat=jerr)                                         
      if(jerr.ne.0) return					                                                 
      allocate(q(1:no,1:nc),stat=jerr)                                          
      if(jerr.ne.0) return					                                                 
      allocate(r(1:no,1:nc),stat=jerr)                                          
      if(jerr.ne.0) return                                                      
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               3967
      exmn=-exmx                                                           3968
      allocate(mm(1:ni),stat=jerr)                                         3969
      if(jerr.ne.0) return                                                 3970
      allocate(is(1:max(nc,ni)),stat=jerr)                                 3971
      if(jerr.ne.0) return                                                 3972
      allocate(sxp(1:no),stat=jerr)                                        3973
      if(jerr.ne.0) return                                                 3974
      allocate(sxpl(1:no),stat=jerr)                                       3975
      if(jerr.ne.0) return                                                 3976
      allocate(ga(1:ni),stat=jerr)                                         3977
      if(jerr.ne.0) return                                                 3978
      allocate(ixx(1:ni),stat=jerr)                                        3979
      if(jerr.ne.0) return                                                 3980
      allocate(gk(1:nc),stat=jerr)                                         3981
      if(jerr.ne.0) return                                                 3982
      allocate(del(1:nc),stat=jerr)                                        3983
      if(jerr.ne.0) return                                                 3984
      allocate(isc(1:nc),stat=jerr)                                        3985
      if(jerr.ne.0) return                                                 3986
      pmax=1.0-pmin                                                        3986
      emin=pmin/pmax                                                       3986
      emax=1.0/emin                                                        3987
      bta=parm                                                             3987
      omb=1.0-bta                                                          3987
      dev1=0.0                                                             3987
      dev0=0.0                                                             3988
21580 do 21581 ic=1,nc                                                     3988
      q0=dot_product(w,y(:,ic))                                            3989
      if(q0 .gt. pmin)goto 21601                                           3989
      jerr =8000+ic                                                        3989
      return                                                               3989
21601 continue                                                             3990
      if(q0 .lt. pmax)goto 21621                                           3990
      jerr =9000+ic                                                        3990
      return                                                               3990
21621 continue                                                             3991
      if(intr .ne. 0)goto 21641                                            3991
      q0=1.0/nc                                                            3991
      b(0,ic)=0.0                                                          3991
      goto 21651                                                           3992
21641 continue                                                             3992
      b(0,ic)=log(q0)                                                      3992
      dev1=dev1-q0*b(0,ic)                                                 3992
21651 continue                                                             3993
21631 continue                                                             3993
      b(1:ni,ic)=0.0                                                       3994
21581 continue                                                             3995
21582 continue                                                             3995
      if(intr.eq.0) dev1=log(float(nc))                                    3995
      ixx=0                                                                3995
      al=0.0                                                               3996
      if(nonzero(no*nc,g) .ne. 0)goto 21671                                3997
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         3997
      sxp=0.0                                                              3998
21680 do 21681 ic=1,nc                                                     3998
      q(:,ic)=exp(b(0,ic))                                                 3998
      sxp=sxp+q(:,ic)                                                      3998
21681 continue                                                             3999
21682 continue                                                             3999
      goto 21691                                                           4000
21671 continue                                                             4000
21700 do 21701 i=1,no                                                      4000
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         4000
21701 continue                                                             4000
21702 continue                                                             4000
      sxp=0.0                                                              4001
      if(intr .ne. 0)goto 21721                                            4001
      b(0,:)=0.0                                                           4001
      goto 21731                                                           4002
21721 continue                                                             4002
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 4002
      if(jerr.ne.0) return                                                 4002
21731 continue                                                             4003
21711 continue                                                             4003
      dev1=0.0                                                             4004
21740 do 21741 ic=1,nc                                                     4004
      q(:,ic)=b(0,ic)+g(:,ic)                                              4005
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             4006
      q(:,ic)=exp(q(:,ic))                                                 4006
      sxp=sxp+q(:,ic)                                                      4007
21741 continue                                                             4008
21742 continue                                                             4008
      sxpl=w*log(sxp)                                                      4008
21750 do 21751 ic=1,nc                                                     4008
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  4008
21751 continue                                                             4009
21752 continue                                                             4009
21691 continue                                                             4010
21661 continue                                                             4010
21760 do 21761 ic=1,nc                                                     4010
21770 do 21771 i=1,no                                                      4010
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               4010
21771 continue                                                             4010
21772 continue                                                             4010
21761 continue                                                             4011
21762 continue                                                             4011
      dev0=dev0+dev1                                                       4013
      alf=1.0                                                              4015
      if(flmin .ge. 1.0)goto 21791                                         4015
      eqs=max(eps,flmin)                                                   4015
      alf=eqs**(1.0/(nlam-1))                                              4015
21791 continue                                                             4016
      m=0                                                                  4016
      mm=0                                                                 4016
      nin=0                                                                4016
      nlp=0                                                                4016
      mnl=min(mnlam,nlam)                                                  4016
      bs=0.0                                                               4016
      shr=shri*dev0                                                        4017
      ga=0.0                                                               4018
21800 do 21801 ic=1,nc                                                     4018
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      4019
21810 do 21811 j=1,ni                                                      4019
      if(ju(j).ne.0) ga(j)=ga(j)+dot_product(r(:,ic),x(:,j))**2            4019
21811 continue                                                             4020
21812 continue                                                             4020
21801 continue                                                             4021
21802 continue                                                             4021
      ga=sqrt(ga)                                                          4022
21820 do 21821 ilm=1,nlam                                                  4022
      al0=al                                                               4023
      if(flmin .lt. 1.0)goto 21841                                         4023
      al=ulam(ilm)                                                         4023
      goto 21831                                                           4024
21841 if(ilm .le. 2)goto 21851                                             4024
      al=al*alf                                                            4024
      goto 21831                                                           4025
21851 if(ilm .ne. 1)goto 21861                                             4025
      al=big                                                               4025
      goto 21871                                                           4026
21861 continue                                                             4026
      al0=0.0                                                              4027
21880 do 21881 j=1,ni                                                      4027
      if(ju(j).eq.0)goto 21881                                             4027
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            4027
21881 continue                                                             4028
21882 continue                                                             4028
      al0=al0/max(bta,1.0d-3)                                              4028
      al=alf*al0                                                           4029
21871 continue                                                             4030
21831 continue                                                             4030
      al2=al*omb                                                           4030
      al1=al*bta                                                           4030
      tlam=bta*(2.0*al-al0)                                                4031
21890 do 21891 k=1,ni                                                      4031
      if(ixx(k).eq.1)goto 21891                                            4031
      if(ju(k).eq.0)goto 21891                                             4032
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     4033
21891 continue                                                             4034
21892 continue                                                             4034
11020 continue                                                             4035
21900 continue                                                             4035
21901 continue                                                             4035
      ix=0                                                                 4035
      jx=ix                                                                4035
      kx=jx                                                                4035
      t=0.0                                                                4036
21910 do 21911 ic=1,nc                                                     4036
      t=max(t,maxval(q(:,ic)*(1.0-q(:,ic)/sxp)/sxp))                       4036
21911 continue                                                             4037
21912 continue                                                             4037
      if(t .ge. eps)goto 21931                                             4037
      kx=1                                                                 4037
      goto 21902                                                           4037
21931 continue                                                             4037
      t=2.0*t                                                              4037
      alt=al1/t                                                            4037
      al2t=al2/t                                                           4038
21940 do 21941 ic=1,nc                                                     4039
      bs(0,ic)=b(0,ic)                                                     4039
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          4040
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)/t                                    4041
      d=0.0                                                                4041
      if(intr.ne.0) d=sum(r(:,ic))                                         4042
      if(d .eq. 0.0)goto 21961                                             4043
      b(0,ic)=b(0,ic)+d                                                    4043
      r(:,ic)=r(:,ic)-d*w                                                  4043
      dlx=max(dlx,d**2)                                                    4044
21961 continue                                                             4045
21941 continue                                                             4046
21942 continue                                                             4046
21970 continue                                                             4046
21971 continue                                                             4046
      nlp=nlp+nc                                                           4046
      dlx=0.0                                                              4047
21980 do 21981 k=1,ni                                                      4047
      if(ixx(k).eq.0)goto 21981                                            4047
      gkn=0.0                                                              4048
21990 do 21991 ic=1,nc                                                     4048
      gk(ic)=dot_product(r(:,ic),x(:,k))+b(k,ic)*xv(k)                     4049
      gkn=gkn+gk(ic)**2                                                    4050
21991 continue                                                             4051
21992 continue                                                             4051
      gkn=sqrt(gkn)                                                        4051
      u=1.0-alt*vp(k)/gkn                                                  4051
      del=b(k,:)                                                           4052
      if(u .gt. 0.0)goto 22011                                             4052
      b(k,:)=0.0                                                           4052
      goto 22021                                                           4053
22011 continue                                                             4053
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                     4054
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),  cl(2,k),vp(k)*al2t,alt*vp(   4056 
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                                 4057
22021 continue                                                             4058
22001 continue                                                             4058
      del=b(k,:)-del                                                       4058
      if(maxval(abs(del)).le.0.0)goto 21981                                4059
22030 do 22031 ic=1,nc                                                     4059
      dlx=max(dlx,xv(k)*del(ic)**2)                                        4060
      r(:,ic)=r(:,ic)-del(ic)*w*x(:,k)                                     4061
22031 continue                                                             4062
22032 continue                                                             4062
      if(mm(k) .ne. 0)goto 22051                                           4062
      nin=nin+1                                                            4063
      if(nin .le. nx)goto 22071                                            4063
      jx=1                                                                 4063
      goto 21982                                                           4063
22071 continue                                                             4064
      mm(k)=nin                                                            4064
      m(nin)=k                                                             4065
22051 continue                                                             4066
21981 continue                                                             4067
21982 continue                                                             4067
      if(jx.gt.0)goto 21972                                                4067
      if(dlx.lt.shr)goto 21972                                             4068
      if(nlp .le. maxit)goto 22091                                         4068
      jerr=-ilm                                                            4068
      return                                                               4068
22091 continue                                                             4069
22100 continue                                                             4069
22101 continue                                                             4069
      nlp=nlp+nc                                                           4069
      dlx=0.0                                                              4070
22110 do 22111 l=1,nin                                                     4070
      k=m(l)                                                               4070
      gkn=0.0                                                              4071
22120 do 22121 ic=1,nc                                                     4071
      gk(ic)=dot_product(r(:,ic),x(:,k))+b(k,ic)*xv(k)                     4072
      gkn=gkn+gk(ic)**2                                                    4073
22121 continue                                                             4074
22122 continue                                                             4074
      gkn=sqrt(gkn)                                                        4074
      u=1.0-alt*vp(k)/gkn                                                  4074
      del=b(k,:)                                                           4075
      if(u .gt. 0.0)goto 22141                                             4075
      b(k,:)=0.0                                                           4075
      goto 22151                                                           4076
22141 continue                                                             4076
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                     4077
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),  cl(2,k),vp(k)*al2t,alt*vp(   4079 
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                                 4080
22151 continue                                                             4081
22131 continue                                                             4081
      del=b(k,:)-del                                                       4081
      if(maxval(abs(del)).le.0.0)goto 22111                                4082
22160 do 22161 ic=1,nc                                                     4082
      dlx=max(dlx,xv(k)*del(ic)**2)                                        4083
      r(:,ic)=r(:,ic)-del(ic)*w*x(:,k)                                     4084
22161 continue                                                             4085
22162 continue                                                             4085
22111 continue                                                             4086
22112 continue                                                             4086
      if(dlx.lt.shr)goto 22102                                             4086
      if(nlp .le. maxit)goto 22181                                         4086
      jerr=-ilm                                                            4086
      return                                                               4086
22181 continue                                                             4088
      goto 22101                                                           4089
22102 continue                                                             4089
      goto 21971                                                           4090
21972 continue                                                             4090
      if(jx.gt.0)goto 21902                                                4091
22190 do 22191 ic=1,nc                                                     4092
      if((b(0,ic)-bs(0,ic))**2.gt.shr) ix=1                                4093
      if(ix .ne. 0)goto 22211                                              4094
22220 do 22221 j=1,nin                                                     4094
      k=m(j)                                                               4095
      if(xv(k)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 22241                   4095
      ix=1                                                                 4095
      goto 22222                                                           4095
22241 continue                                                             4097
22221 continue                                                             4098
22222 continue                                                             4098
22211 continue                                                             4099
22250 do 22251 i=1,no                                                      4099
      fi=b(0,ic)+g(i,ic)                                                   4101
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)))         4102
      fi=min(max(exmn,fi),exmx)                                            4102
      sxp(i)=sxp(i)-q(i,ic)                                                4103
      q(i,ic)=min(max(emin*sxp(i),exp(fi)),emax*sxp(i))                    4104
      sxp(i)=sxp(i)+q(i,ic)                                                4105
22251 continue                                                             4106
22252 continue                                                             4106
22191 continue                                                             4107
22192 continue                                                             4107
      s=-sum(b(0,:))/nc                                                    4107
      b(0,:)=b(0,:)+s                                                      4108
      if(jx.gt.0)goto 21902                                                4109
      if(ix .ne. 0)goto 22271                                              4110
22280 do 22281 k=1,ni                                                      4110
      if(ixx(k).eq.1)goto 22281                                            4110
      if(ju(k).eq.0)goto 22281                                             4110
      ga(k)=0.0                                                            4110
22281 continue                                                             4111
22282 continue                                                             4111
22290 do 22291 ic=1,nc                                                     4111
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      4112
22300 do 22301 k=1,ni                                                      4112
      if(ixx(k).eq.1)goto 22301                                            4112
      if(ju(k).eq.0)goto 22301                                             4113
      ga(k)=ga(k)+dot_product(r(:,ic),x(:,k))**2                           4114
22301 continue                                                             4115
22302 continue                                                             4115
22291 continue                                                             4116
22292 continue                                                             4116
      ga=sqrt(ga)                                                          4117
22310 do 22311 k=1,ni                                                      4117
      if(ixx(k).eq.1)goto 22311                                            4117
      if(ju(k).eq.0)goto 22311                                             4118
      if(ga(k) .le. al1*vp(k))goto 22331                                   4118
      ixx(k)=1                                                             4118
      ix=1                                                                 4118
22331 continue                                                             4119
22311 continue                                                             4120
22312 continue                                                             4120
      if(ix.eq.1) go to 11020                                              4121
      goto 21902                                                           4122
22271 continue                                                             4123
      goto 21901                                                           4124
21902 continue                                                             4124
      if(kx .le. 0)goto 22351                                              4124
      jerr=-20000-ilm                                                      4124
      goto 21822                                                           4124
22351 continue                                                             4125
      if(jx .le. 0)goto 22371                                              4125
      jerr=-10000-ilm                                                      4125
      goto 21822                                                           4125
22371 continue                                                             4125
      devi=0.0                                                             4126
22380 do 22381 ic=1,nc                                                     4127
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          4127
      a0(ic,ilm)=b(0,ic)                                                   4128
22390 do 22391 i=1,no                                                      4128
      if(y(i,ic).le.0.0)goto 22391                                         4129
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           4130
22391 continue                                                             4131
22392 continue                                                             4131
22381 continue                                                             4132
22382 continue                                                             4132
      kin(ilm)=nin                                                         4132
      alm(ilm)=al                                                          4132
      lmu=ilm                                                              4133
      dev(ilm)=(dev1-devi)/dev0                                            4134
      if(ilm.lt.mnl)goto 21821                                             4134
      if(flmin.ge.1.0)goto 21821                                           4135
      me=0                                                                 4135
22400 do 22401 j=1,nin                                                     4135
      if(a(j,1,ilm).ne.0.0) me=me+1                                        4135
22401 continue                                                             4135
22402 continue                                                             4135
      if(me.gt.ne)goto 21822                                               4136
      if(dev(ilm).gt.devmax)goto 21822                                     4136
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 21822                             4137
21821 continue                                                             4138
21822 continue                                                             4138
      g=log(q)                                                             4138
22410 do 22411 i=1,no                                                      4138
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         4138
22411 continue                                                             4139
22412 continue                                                             4139
      deallocate(sxp,b,bs,r,q,mm,is,ga,ixx,gk,del,sxpl)                    4140
      return                                                               4141
      end                                                                  4142
      subroutine multsprlognetn(parm,no,ni,nc,x,ix,jx,y,g,w,ju,vp,cl,ne,   4144 
     *nx,nlam,  flmin,ulam,shri,intr,maxit,xv,xb,xs,lmu,a0,a,m,kin,dev0,
     *dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   4145
      double precision x(*),y(no,nc),g(no,nc),w(no),vp(ni)                 4146
      double precision ulam(nlam),xb(ni),xs(ni),xv(ni)                     4147
      double precision a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(   4148 
     *2,ni)
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           4149
      double precision, dimension (:,:), allocatable :: q,r,b,bs                
      double precision, dimension (:), allocatable :: sxp,sxpl,ga,gk            
      double precision, dimension (:), allocatable :: del,sc,svr                
      integer, dimension (:), allocatable :: mm,is,iy,isc                       
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      if(jerr.ne.0) return                                                      
      allocate(bs(0:ni,1:nc),stat=jerr)                                         
      if(jerr.ne.0) return					                                                 
      allocate(q(1:no,1:nc),stat=jerr)                                          
      if(jerr.ne.0) return					                                                 
      allocate(r(1:no,1:nc),stat=jerr)                                          
      if(jerr.ne.0) return					                                                 
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               4163
      exmn=-exmx                                                           4164
      allocate(mm(1:ni),stat=jerr)                                         4165
      if(jerr.ne.0) return                                                 4166
      allocate(ga(1:ni),stat=jerr)                                         4167
      if(jerr.ne.0) return                                                 4168
      allocate(gk(1:nc),stat=jerr)                                         4169
      if(jerr.ne.0) return                                                 4170
      allocate(del(1:nc),stat=jerr)                                        4171
      if(jerr.ne.0) return                                                 4172
      allocate(iy(1:ni),stat=jerr)                                         4173
      if(jerr.ne.0) return                                                 4174
      allocate(is(1:max(nc,ni)),stat=jerr)                                 4175
      if(jerr.ne.0) return                                                 4176
      allocate(sxp(1:no),stat=jerr)                                        4177
      if(jerr.ne.0) return                                                 4178
      allocate(sxpl(1:no),stat=jerr)                                       4179
      if(jerr.ne.0) return                                                 4180
      allocate(svr(1:nc),stat=jerr)                                        4181
      if(jerr.ne.0) return                                                 4182
      allocate(sc(1:no),stat=jerr)                                         4183
      if(jerr.ne.0) return                                                 4184
      allocate(isc(1:nc),stat=jerr)                                        4185
      if(jerr.ne.0) return                                                 4186
      pmax=1.0-pmin                                                        4186
      emin=pmin/pmax                                                       4186
      emax=1.0/emin                                                        4187
      bta=parm                                                             4187
      omb=1.0-bta                                                          4187
      dev1=0.0                                                             4187
      dev0=0.0                                                             4188
22420 do 22421 ic=1,nc                                                     4188
      q0=dot_product(w,y(:,ic))                                            4189
      if(q0 .gt. pmin)goto 22441                                           4189
      jerr =8000+ic                                                        4189
      return                                                               4189
22441 continue                                                             4190
      if(q0 .lt. pmax)goto 22461                                           4190
      jerr =9000+ic                                                        4190
      return                                                               4190
22461 continue                                                             4191
      b(1:ni,ic)=0.0                                                       4192
      if(intr .ne. 0)goto 22481                                            4192
      q0=1.0/nc                                                            4192
      b(0,ic)=0.0                                                          4192
      goto 22491                                                           4193
22481 continue                                                             4193
      b(0,ic)=log(q0)                                                      4193
      dev1=dev1-q0*b(0,ic)                                                 4193
22491 continue                                                             4194
22471 continue                                                             4194
22421 continue                                                             4195
22422 continue                                                             4195
      if(intr.eq.0) dev1=log(float(nc))                                    4195
      iy=0                                                                 4195
      al=0.0                                                               4196
      if(nonzero(no*nc,g) .ne. 0)goto 22511                                4197
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         4197
      sxp=0.0                                                              4198
22520 do 22521 ic=1,nc                                                     4198
      q(:,ic)=exp(b(0,ic))                                                 4198
      sxp=sxp+q(:,ic)                                                      4198
22521 continue                                                             4199
22522 continue                                                             4199
      goto 22531                                                           4200
22511 continue                                                             4200
22540 do 22541 i=1,no                                                      4200
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         4200
22541 continue                                                             4200
22542 continue                                                             4200
      sxp=0.0                                                              4201
      if(intr .ne. 0)goto 22561                                            4201
      b(0,:)=0.0                                                           4201
      goto 22571                                                           4202
22561 continue                                                             4202
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 4202
      if(jerr.ne.0) return                                                 4202
22571 continue                                                             4203
22551 continue                                                             4203
      dev1=0.0                                                             4204
22580 do 22581 ic=1,nc                                                     4204
      q(:,ic)=b(0,ic)+g(:,ic)                                              4205
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             4206
      q(:,ic)=exp(q(:,ic))                                                 4206
      sxp=sxp+q(:,ic)                                                      4207
22581 continue                                                             4208
22582 continue                                                             4208
      sxpl=w*log(sxp)                                                      4208
22590 do 22591 ic=1,nc                                                     4208
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  4208
22591 continue                                                             4209
22592 continue                                                             4209
22531 continue                                                             4210
22501 continue                                                             4210
22600 do 22601 ic=1,nc                                                     4210
22610 do 22611 i=1,no                                                      4210
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               4210
22611 continue                                                             4210
22612 continue                                                             4210
22601 continue                                                             4211
22602 continue                                                             4211
      dev0=dev0+dev1                                                       4213
      alf=1.0                                                              4215
      if(flmin .ge. 1.0)goto 22631                                         4215
      eqs=max(eps,flmin)                                                   4215
      alf=eqs**(1.0/(nlam-1))                                              4215
22631 continue                                                             4216
      m=0                                                                  4216
      mm=0                                                                 4216
      nin=0                                                                4216
      nlp=0                                                                4216
      mnl=min(mnlam,nlam)                                                  4216
      bs=0.0                                                               4217
      shr=shri*dev0                                                        4217
      ga=0.0                                                               4218
22640 do 22641 ic=1,nc                                                     4218
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      4218
      svr(ic)=sum(r(:,ic))                                                 4219
22650 do 22651 j=1,ni                                                      4219
      if(ju(j).eq.0)goto 22651                                             4220
      jb=ix(j)                                                             4220
      je=ix(j+1)-1                                                         4221
      gj=dot_product(r(jx(jb:je),ic),x(jb:je))                             4222
      ga(j)=ga(j)+((gj-svr(ic)*xb(j))/xs(j))**2                            4223
22651 continue                                                             4224
22652 continue                                                             4224
22641 continue                                                             4225
22642 continue                                                             4225
      ga=sqrt(ga)                                                          4226
22660 do 22661 ilm=1,nlam                                                  4226
      al0=al                                                               4227
      if(flmin .lt. 1.0)goto 22681                                         4227
      al=ulam(ilm)                                                         4227
      goto 22671                                                           4228
22681 if(ilm .le. 2)goto 22691                                             4228
      al=al*alf                                                            4228
      goto 22671                                                           4229
22691 if(ilm .ne. 1)goto 22701                                             4229
      al=big                                                               4229
      goto 22711                                                           4230
22701 continue                                                             4230
      al0=0.0                                                              4231
22720 do 22721 j=1,ni                                                      4231
      if(ju(j).eq.0)goto 22721                                             4231
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            4231
22721 continue                                                             4232
22722 continue                                                             4232
      al0=al0/max(bta,1.0d-3)                                              4232
      al=alf*al0                                                           4233
22711 continue                                                             4234
22671 continue                                                             4234
      al2=al*omb                                                           4234
      al1=al*bta                                                           4234
      tlam=bta*(2.0*al-al0)                                                4235
22730 do 22731 k=1,ni                                                      4235
      if(iy(k).eq.1)goto 22731                                             4235
      if(ju(k).eq.0)goto 22731                                             4236
      if(ga(k).gt.tlam*vp(k)) iy(k)=1                                      4237
22731 continue                                                             4238
22732 continue                                                             4238
11020 continue                                                             4239
22740 continue                                                             4239
22741 continue                                                             4239
      ixx=0                                                                4239
      jxx=ixx                                                              4239
      kxx=jxx                                                              4239
      t=0.0                                                                4240
22750 do 22751 ic=1,nc                                                     4240
      t=max(t,maxval(q(:,ic)*(1.0-q(:,ic)/sxp)/sxp))                       4240
22751 continue                                                             4241
22752 continue                                                             4241
      if(t .ge. eps)goto 22771                                             4241
      kxx=1                                                                4241
      goto 22742                                                           4241
22771 continue                                                             4241
      t=2.0*t                                                              4241
      alt=al1/t                                                            4241
      al2t=al2/t                                                           4242
22780 do 22781 ic=1,nc                                                     4242
      bs(0,ic)=b(0,ic)                                                     4242
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          4243
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)/t                                    4243
      svr(ic)=sum(r(:,ic))                                                 4244
      if(intr .eq. 0)goto 22801                                            4244
      b(0,ic)=b(0,ic)+svr(ic)                                              4244
      r(:,ic)=r(:,ic)-svr(ic)*w                                            4245
      dlx=max(dlx,svr(ic)**2)                                              4246
22801 continue                                                             4247
22781 continue                                                             4248
22782 continue                                                             4248
22810 continue                                                             4248
22811 continue                                                             4248
      nlp=nlp+nc                                                           4248
      dlx=0.0                                                              4249
22820 do 22821 k=1,ni                                                      4249
      if(iy(k).eq.0)goto 22821                                             4250
      jb=ix(k)                                                             4250
      je=ix(k+1)-1                                                         4250
      del=b(k,:)                                                           4250
      gkn=0.0                                                              4251
22830 do 22831 ic=1,nc                                                     4252
      u=(dot_product(r(jx(jb:je),ic),x(jb:je))-svr(ic)*xb(k))/xs(k)        4253
      gk(ic)=u+del(ic)*xv(k)                                               4253
      gkn=gkn+gk(ic)**2                                                    4254
22831 continue                                                             4255
22832 continue                                                             4255
      gkn=sqrt(gkn)                                                        4255
      u=1.0-alt*vp(k)/gkn                                                  4256
      if(u .gt. 0.0)goto 22851                                             4256
      b(k,:)=0.0                                                           4256
      goto 22861                                                           4257
22851 continue                                                             4258
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                     4259
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),cl(2,k),  vp(k)*al2t,alt*vp(   4261 
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                                 4262
22861 continue                                                             4263
22841 continue                                                             4263
      del=b(k,:)-del                                                       4263
      if(maxval(abs(del)).le.0.0)goto 22821                                4264
22870 do 22871 ic=1,nc                                                     4264
      dlx=max(dlx,xv(k)*del(ic)**2)                                        4265
      r(jx(jb:je),ic)=r(jx(jb:je),ic)  -del(ic)*w(jx(jb:je))*(x(jb:je)-x   4267 
     *b(k))/xs(k)
22871 continue                                                             4268
22872 continue                                                             4268
      if(mm(k) .ne. 0)goto 22891                                           4268
      nin=nin+1                                                            4269
      if(nin .le. nx)goto 22911                                            4269
      jxx=1                                                                4269
      goto 22822                                                           4269
22911 continue                                                             4270
      mm(k)=nin                                                            4270
      m(nin)=k                                                             4271
22891 continue                                                             4272
22821 continue                                                             4273
22822 continue                                                             4273
      if(jxx.gt.0)goto 22812                                               4274
      if(dlx.lt.shr)goto 22812                                             4274
      if(nlp .le. maxit)goto 22931                                         4274
      jerr=-ilm                                                            4274
      return                                                               4274
22931 continue                                                             4275
22940 continue                                                             4275
22941 continue                                                             4275
      nlp=nlp+nc                                                           4275
      dlx=0.0                                                              4276
22950 do 22951 l=1,nin                                                     4276
      k=m(l)                                                               4276
      jb=ix(k)                                                             4276
      je=ix(k+1)-1                                                         4276
      del=b(k,:)                                                           4276
      gkn=0.0                                                              4277
22960 do 22961 ic=1,nc                                                     4278
      u=(dot_product(r(jx(jb:je),ic),x(jb:je))  -svr(ic)*xb(k))/xs(k)      4280
      gk(ic)=u+del(ic)*xv(k)                                               4280
      gkn=gkn+gk(ic)**2                                                    4281
22961 continue                                                             4282
22962 continue                                                             4282
      gkn=sqrt(gkn)                                                        4282
      u=1.0-alt*vp(k)/gkn                                                  4283
      if(u .gt. 0.0)goto 22981                                             4283
      b(k,:)=0.0                                                           4283
      goto 22991                                                           4284
22981 continue                                                             4285
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                     4286
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),cl(2,k),  vp(k)*al2t,alt*vp(   4288 
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                                 4289
22991 continue                                                             4290
22971 continue                                                             4290
      del=b(k,:)-del                                                       4290
      if(maxval(abs(del)).le.0.0)goto 22951                                4291
23000 do 23001 ic=1,nc                                                     4291
      dlx=max(dlx,xv(k)*del(ic)**2)                                        4292
      r(jx(jb:je),ic)=r(jx(jb:je),ic)  -del(ic)*w(jx(jb:je))*(x(jb:je)-x   4294 
     *b(k))/xs(k)
23001 continue                                                             4295
23002 continue                                                             4295
22951 continue                                                             4296
22952 continue                                                             4296
      if(dlx.lt.shr)goto 22942                                             4296
      if(nlp .le. maxit)goto 23021                                         4296
      jerr=-ilm                                                            4296
      return                                                               4296
23021 continue                                                             4298
      goto 22941                                                           4299
22942 continue                                                             4299
      goto 22811                                                           4300
22812 continue                                                             4300
      if(jxx.gt.0)goto 22742                                               4301
23030 do 23031 ic=1,nc                                                     4302
      if((b(0,ic)-bs(0,ic))**2.gt.shr) ixx=1                               4303
      if(ixx .ne. 0)goto 23051                                             4304
23060 do 23061 j=1,nin                                                     4304
      k=m(j)                                                               4305
      if(xv(k)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 23081                   4305
      ixx=1                                                                4305
      goto 23062                                                           4305
23081 continue                                                             4307
23061 continue                                                             4308
23062 continue                                                             4308
23051 continue                                                             4309
      sc=b(0,ic)+g(:,ic)                                                   4309
      b0=0.0                                                               4310
23090 do 23091 j=1,nin                                                     4310
      l=m(j)                                                               4310
      jb=ix(l)                                                             4310
      je=ix(l+1)-1                                                         4311
      sc(jx(jb:je))=sc(jx(jb:je))+b(l,ic)*x(jb:je)/xs(l)                   4312
      b0=b0-b(l,ic)*xb(l)/xs(l)                                            4313
23091 continue                                                             4314
23092 continue                                                             4314
      sc=min(max(exmn,sc+b0),exmx)                                         4315
      sxp=sxp-q(:,ic)                                                      4316
      q(:,ic)=min(max(emin*sxp,exp(sc)),emax*sxp)                          4317
      sxp=sxp+q(:,ic)                                                      4318
23031 continue                                                             4319
23032 continue                                                             4319
      s=sum(b(0,:))/nc                                                     4319
      b(0,:)=b(0,:)-s                                                      4320
      if(jxx.gt.0)goto 22742                                               4321
      if(ixx .ne. 0)goto 23111                                             4322
23120 do 23121 j=1,ni                                                      4322
      if(iy(j).eq.1)goto 23121                                             4322
      if(ju(j).eq.0)goto 23121                                             4322
      ga(j)=0.0                                                            4322
23121 continue                                                             4323
23122 continue                                                             4323
23130 do 23131 ic=1,nc                                                     4323
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      4324
23140 do 23141 j=1,ni                                                      4324
      if(iy(j).eq.1)goto 23141                                             4324
      if(ju(j).eq.0)goto 23141                                             4325
      jb=ix(j)                                                             4325
      je=ix(j+1)-1                                                         4326
      gj=dot_product(r(jx(jb:je),ic),x(jb:je))                             4327
      ga(j)=ga(j)+((gj-svr(ic)*xb(j))/xs(j))**2                            4328
23141 continue                                                             4329
23142 continue                                                             4329
23131 continue                                                             4330
23132 continue                                                             4330
      ga=sqrt(ga)                                                          4331
23150 do 23151 k=1,ni                                                      4331
      if(iy(k).eq.1)goto 23151                                             4331
      if(ju(k).eq.0)goto 23151                                             4332
      if(ga(k) .le. al1*vp(k))goto 23171                                   4332
      iy(k)=1                                                              4332
      ixx=1                                                                4332
23171 continue                                                             4333
23151 continue                                                             4334
23152 continue                                                             4334
      if(ixx.eq.1) go to 11020                                             4335
      goto 22742                                                           4336
23111 continue                                                             4337
      goto 22741                                                           4338
22742 continue                                                             4338
      if(kxx .le. 0)goto 23191                                             4338
      jerr=-20000-ilm                                                      4338
      goto 22662                                                           4338
23191 continue                                                             4339
      if(jxx .le. 0)goto 23211                                             4339
      jerr=-10000-ilm                                                      4339
      goto 22662                                                           4339
23211 continue                                                             4339
      devi=0.0                                                             4340
23220 do 23221 ic=1,nc                                                     4341
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          4341
      a0(ic,ilm)=b(0,ic)                                                   4342
23230 do 23231 i=1,no                                                      4342
      if(y(i,ic).le.0.0)goto 23231                                         4343
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           4344
23231 continue                                                             4345
23232 continue                                                             4345
23221 continue                                                             4346
23222 continue                                                             4346
      kin(ilm)=nin                                                         4346
      alm(ilm)=al                                                          4346
      lmu=ilm                                                              4347
      dev(ilm)=(dev1-devi)/dev0                                            4348
      if(ilm.lt.mnl)goto 22661                                             4348
      if(flmin.ge.1.0)goto 22661                                           4349
      me=0                                                                 4349
23240 do 23241 j=1,nin                                                     4349
      if(a(j,1,ilm).ne.0.0) me=me+1                                        4349
23241 continue                                                             4349
23242 continue                                                             4349
      if(me.gt.ne)goto 22662                                               4350
      if(dev(ilm).gt.devmax)goto 22662                                     4350
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 22662                             4351
22661 continue                                                             4352
22662 continue                                                             4352
      g=log(q)                                                             4352
23250 do 23251 i=1,no                                                      4352
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         4352
23251 continue                                                             4353
23252 continue                                                             4353
      deallocate(sxp,b,bs,r,q,mm,is,sc,ga,iy,gk,del,sxpl)                  4354
      return                                                               4355
      end                                                                  4356
      subroutine psort7 (v,a,ii,jj)                                             
      implicit double precision(a-h,o-z)                                        
c                                                                               
c     puts into a the permutation vector which sorts v into                     
c     increasing order. the array v is not modified.                            
c     only elements from ii to jj are considered.                               
c     arrays iu(k) and il(k) permit sorting up to 2**(k+1)-1 elements           
c                                                                               
c     this is a modification of cacm algorithm #347 by r. c. singleton,         
c     which is a modified hoare quicksort.                                      
c                                                                               
      dimension a(jj),v(jj),iu(20),il(20)                                       
      integer t,tt                                                              
      integer a                                                                 
      double precision v                                                        
      m=1                                                                       
      i=ii                                                                      
      j=jj                                                                      
 10   if (i.ge.j) go to 80                                                      
 20   k=i                                                                       
      ij=(j+i)/2                                                                
      t=a(ij)                                                                   
      vt=v(t)                                                                   
      if (v(a(i)).le.vt) go to 30                                               
      a(ij)=a(i)                                                                
      a(i)=t                                                                    
      t=a(ij)                                                                   
      vt=v(t)                                                                   
 30   l=j                                                                       
      if (v(a(j)).ge.vt) go to 50                                               
      a(ij)=a(j)                                                                
      a(j)=t                                                                    
      t=a(ij)                                                                   
      vt=v(t)                                                                   
      if (v(a(i)).le.vt) go to 50                                               
      a(ij)=a(i)                                                                
      a(i)=t                                                                    
      t=a(ij)                                                                   
      vt=v(t)                                                                   
      go to 50                                                                  
 40   a(l)=a(k)                                                                 
      a(k)=tt                                                                   
 50   l=l-1                                                                     
      if (v(a(l)).gt.vt) go to 50                                               
      tt=a(l)                                                                   
      vtt=v(tt)                                                                 
 60   k=k+1                                                                     
      if (v(a(k)).lt.vt) go to 60                                               
      if (k.le.l) go to 40                                                      
      if (l-i.le.j-k) go to 70                                                  
      il(m)=i                                                                   
      iu(m)=l                                                                   
      i=k                                                                       
      m=m+1                                                                     
      go to 90                                                                  
 70   il(m)=k                                                                   
      iu(m)=j                                                                   
      j=l                                                                       
      m=m+1                                                                     
      go to 90                                                                  
 80   m=m-1                                                                     
      if (m.eq.0) return                                                        
      i=il(m)                                                                   
      j=iu(m)                                                                   
 90   if (j-i.gt.10) go to 20                                                   
      if (i.eq.ii) go to 10                                                     
      i=i-1                                                                     
 100  i=i+1                                                                     
      if (i.eq.j) go to 80                                                      
      t=a(i+1)                                                                  
      vt=v(t)                                                                   
      if (v(a(i)).le.vt) go to 100                                              
      k=i                                                                       
 110  a(k+1)=a(k)                                                               
      k=k-1                                                                     
      if (vt.lt.v(a(k))) go to 110                                              
      a(k+1)=t                                                                  
      go to 100                                                                 
      end                                                                       
