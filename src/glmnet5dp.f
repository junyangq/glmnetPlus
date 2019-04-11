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
     *lp,jerr,res)
      implicit double precision(a-h,o-z)                                    790
      double precision x(no,ni),y(no),w(no),vp(ni),ca(nx,nlam),cl(2,ni)     791
      double precision beta0(ni),res(no,nlam)                               792
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
     *sd,intr,maxit,  beta0,isg,plam,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr,r
     *es)
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
     *jerr,res)
      implicit double precision(a-h,o-z)                                    978
      double precision vp(ni),x(no,ni),y(no),w(no),ulam(nlam),cl(2,ni)      979
      double precision beta0(ni),res(no,nlam)                               980
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
     *t,  beta0,isg,plam,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr,res)
      if(jerr.gt.0) return                                                 1008
10640 do 10641 k=1,lmu                                                     1008
      alm(k)=ys*alm(k)                                                     1008
      nk=nin(k)                                                            1008
      res(:,k)=ys*sqrt(float(no))*res(:,k)                                 1009
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
     *r,maxit,  beta0,isg,plam,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr,res)
      implicit double precision(a-h,o-z)                                   1049
      double precision vp(ni),y(no),x(no,ni),ulam(nlam),ao(nx,nlam)        1050
      double precision rsqo(nlam),almo(nlam),xv(ni)                        1051
      double precision cl(2,ni),beta0(ni),res(no,nlam)                     1052
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
      if(abs(beta0(k)) .le. 1e-24)goto 10951                               1085
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
      res(:,m)=y                                                           1140
      if(m.lt.mnl)goto 10841                                               1140
      if(flmin.ge.1.0)goto 10841                                           1141
      me=0                                                                 1141
11190 do 11191 j=1,nin                                                     1141
      if(ao(j,m).ne.0.0) me=me+1                                           1141
11191 continue                                                             1141
11192 continue                                                             1141
      if(me.gt.ne)goto 10842                                               1142
      if(rsq-rsq0.lt.sml*rsq)goto 10842                                    1142
      if(rsq.gt.rsqmax)goto 10842                                          1143
10841 continue                                                             1144
10842 continue                                                             1144
      deallocate(a,mm,g,ix)                                                1145
      return                                                               1146
      end                                                                  1147
      subroutine chkvars(no,ni,x,ju)                                       1148
      implicit double precision(a-h,o-z)                                   1149
      double precision x(no,ni)                                            1149
      integer ju(ni)                                                       1150
11200 do 11201 j=1,ni                                                      1150
      ju(j)=0                                                              1150
      t=x(1,j)                                                             1151
11210 do 11211 i=2,no                                                      1151
      if(x(i,j).eq.t)goto 11211                                            1151
      ju(j)=1                                                              1151
      goto 11212                                                           1151
11211 continue                                                             1152
11212 continue                                                             1152
11201 continue                                                             1153
11202 continue                                                             1153
      return                                                               1154
      end                                                                  1155
      subroutine uncomp(ni,ca,ia,nin,a)                                    1156
      implicit double precision(a-h,o-z)                                   1157
      double precision ca(*),a(ni)                                         1157
      integer ia(*)                                                        1158
      a=0.0                                                                1158
      if(nin.gt.0) a(ia(1:nin))=ca(1:nin)                                  1159
      return                                                               1160
      end                                                                  1161
      subroutine modval(a0,ca,ia,nin,n,x,f)                                1162
      implicit double precision(a-h,o-z)                                   1163
      double precision ca(nin),x(n,*),f(n)                                 1163
      integer ia(nin)                                                      1164
      f=a0                                                                 1164
      if(nin.le.0) return                                                  1165
11220 do 11221 i=1,n                                                       1165
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                      1165
11221 continue                                                             1166
11222 continue                                                             1166
      return                                                               1167
      end                                                                  1168
      subroutine spelnet  (ka,parm,no,ni,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam   1171 
     *,flmin,ulam,thr,isd,intr,  maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr
     *)
      implicit double precision(a-h,o-z)                                   1172
      double precision x(*),y(no),w(no),vp(ni),ulam(nlam),cl(2,ni)         1173
      double precision ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)            1174
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1175
      double precision, dimension (:), allocatable :: vq;                       
      if(maxval(vp) .gt. 0.0)goto 11241                                    1178
      jerr=10000                                                           1178
      return                                                               1178
11241 continue                                                             1179
      allocate(vq(1:ni),stat=jerr)                                         1179
      if(jerr.ne.0) return                                                 1180
      vq=max(0d0,vp)                                                       1180
      vq=vq*ni/sum(vq)                                                     1181
      if(ka .ne. 1)goto 11261                                              1182
      call spelnetu  (parm,no,ni,x,ix,jx,y,w,jd,vq,cl,ne,nx,nlam,flmin,u   1185 
     *lam,thr,isd,  intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 11271                                                           1186
11261 continue                                                             1187
      call spelnetn (parm,no,ni,x,ix,jx,y,w,jd,vq,cl,ne,nx,nlam,flmin,ul   1190 
     *am,thr,isd,intr,  maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
11271 continue                                                             1191
11251 continue                                                             1191
      deallocate(vq)                                                       1192
      return                                                               1193
      end                                                                  1194
      subroutine spelnetu  (parm,no,ni,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam,f   1197 
     *lmin,ulam,thr,isd,intr,  maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   1198
      double precision x(*),y(no),w(no),vp(ni),ulam(nlam),cl(2,ni)         1199
      double precision ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)            1200
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1201
      double precision, dimension (:), allocatable :: xm,xs,g,xv,vlam           
      integer, dimension (:), allocatable :: ju                                 
      allocate(g(1:ni),stat=jerr)                                          1206
      if(jerr.ne.0) return                                                 1207
      allocate(xm(1:ni),stat=jerr)                                         1208
      if(jerr.ne.0) return                                                 1209
      allocate(xs(1:ni),stat=jerr)                                         1210
      if(jerr.ne.0) return                                                 1211
      allocate(ju(1:ni),stat=jerr)                                         1212
      if(jerr.ne.0) return                                                 1213
      allocate(xv(1:ni),stat=jerr)                                         1214
      if(jerr.ne.0) return                                                 1215
      allocate(vlam(1:nlam),stat=jerr)                                     1216
      if(jerr.ne.0) return                                                 1217
      call spchkvars(no,ni,x,ix,ju)                                        1218
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1219
      if(maxval(ju) .gt. 0)goto 11291                                      1219
      jerr=7777                                                            1219
      return                                                               1219
11291 continue                                                             1220
      call spstandard(no,ni,x,ix,jx,y,w,ju,isd,intr,g,xm,xs,ym,ys,xv,jer   1221 
     *r)
      if(jerr.ne.0) return                                                 1222
      cl=cl/ys                                                             1222
      if(isd .le. 0)goto 11311                                             1222
11320 do 11321 j=1,ni                                                      1222
      cl(:,j)=cl(:,j)*xs(j)                                                1222
11321 continue                                                             1222
11322 continue                                                             1222
11311 continue                                                             1223
      if(flmin.ge.1.0) vlam=ulam/ys                                        1224
      call spelnet1(parm,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,cl,nlam,flmin,vla   1226 
     *m,thr,maxit,  xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 1227
11330 do 11331 k=1,lmu                                                     1227
      alm(k)=ys*alm(k)                                                     1227
      nk=nin(k)                                                            1228
11340 do 11341 l=1,nk                                                      1228
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                         1228
11341 continue                                                             1228
11342 continue                                                             1228
      a0(k)=0.0                                                            1229
      if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))          1230
11331 continue                                                             1231
11332 continue                                                             1231
      deallocate(xm,xs,g,ju,xv,vlam)                                       1232
      return                                                               1233
      end                                                                  1234
      subroutine spstandard (no,ni,x,ix,jx,y,w,ju,isd,intr,g,xm,xs,ym,ys   1235 
     *,xv,jerr)
      implicit double precision(a-h,o-z)                                   1236
      double precision x(*),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)         1237
      integer ix(*),jx(*),ju(ni)                                           1238
      w=w/sum(w)                                                           1239
      if(intr .ne. 0)goto 11361                                            1239
      ym=0.0                                                               1240
      ys=sqrt(dot_product(w,y**2)-dot_product(w,y)**2)                     1240
      y=y/ys                                                               1241
11370 do 11371 j=1,ni                                                      1241
      if(ju(j).eq.0)goto 11371                                             1241
      xm(j)=0.0                                                            1241
      jb=ix(j)                                                             1241
      je=ix(j+1)-1                                                         1242
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)                          1243
      if(isd .eq. 0)goto 11391                                             1243
      xbq=dot_product(w(jx(jb:je)),x(jb:je))**2                            1243
      vc=xv(j)-xbq                                                         1244
      xs(j)=sqrt(vc)                                                       1244
      xv(j)=1.0+xbq/vc                                                     1245
      goto 11401                                                           1246
11391 continue                                                             1246
      xs(j)=1.0                                                            1246
11401 continue                                                             1247
11381 continue                                                             1247
11371 continue                                                             1248
11372 continue                                                             1248
      goto 11411                                                           1249
11361 continue                                                             1250
11420 do 11421 j=1,ni                                                      1250
      if(ju(j).eq.0)goto 11421                                             1251
      jb=ix(j)                                                             1251
      je=ix(j+1)-1                                                         1251
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1252
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 1253
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       1254
11421 continue                                                             1255
11422 continue                                                             1255
      if(isd .ne. 0)goto 11441                                             1255
      xs=1.0                                                               1255
      goto 11451                                                           1255
11441 continue                                                             1255
      xv=1.0                                                               1255
11451 continue                                                             1256
11431 continue                                                             1256
      ym=dot_product(w,y)                                                  1256
      y=y-ym                                                               1256
      ys=sqrt(dot_product(w,y**2))                                         1256
      y=y/ys                                                               1257
11411 continue                                                             1258
11351 continue                                                             1258
      g=0.0                                                                1259
11460 do 11461 j=1,ni                                                      1259
      if(ju(j).eq.0)goto 11461                                             1259
      jb=ix(j)                                                             1259
      je=ix(j+1)-1                                                         1260
      g(j)=dot_product(w(jx(jb:je))*y(jx(jb:je)),x(jb:je))/xs(j)           1261
11461 continue                                                             1262
11462 continue                                                             1262
      return                                                               1263
      end                                                                  1264
      subroutine spelnet1(beta,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,cl,nlam,flm   1266 
     *in,ulam,  thr,maxit,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      implicit double precision(a-h,o-z)                                   1267
      double precision g(ni),vp(ni),x(*),ulam(nlam),w(no)                  1268
      double precision ao(nx,nlam),rsqo(nlam),almo(nlam)                   1269
      double precision xm(ni),xs(ni),xv(ni),cl(2,ni)                       1270
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                          1271
      double precision, dimension (:), allocatable :: a,da                      
      integer, dimension (:), allocatable :: mm                                 
      double precision, dimension (:,:), allocatable :: c                       
      allocate(c(1:ni,1:nx),stat=jerr)                                          
      if(jerr.ne.0) return;                                                     
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)               1279
      allocate(a(1:ni),stat=jerr)                                          1280
      if(jerr.ne.0) return                                                 1281
      allocate(mm(1:ni),stat=jerr)                                         1282
      if(jerr.ne.0) return                                                 1283
      allocate(da(1:ni),stat=jerr)                                         1284
      if(jerr.ne.0) return                                                 1285
      bta=beta                                                             1285
      omb=1.0-bta                                                          1287
      alm=0.0                                                              1287
      alf=1.0                                                              1289
      if(flmin .ge. 1.0)goto 11481                                         1289
      eqs=max(eps,flmin)                                                   1289
      alf=eqs**(1.0/(nlam-1))                                              1289
11481 continue                                                             1290
      rsq=0.0                                                              1290
      a=0.0                                                                1290
      mm=0                                                                 1290
      nlp=0                                                                1290
      nin=nlp                                                              1290
      iz=0                                                                 1290
      mnl=min(mnlam,nlam)                                                  1291
11490 do 11491 m=1,nlam                                                    1292
      if(flmin .lt. 1.0)goto 11511                                         1292
      alm=ulam(m)                                                          1292
      goto 11501                                                           1293
11511 if(m .le. 2)goto 11521                                               1293
      alm=alm*alf                                                          1293
      goto 11501                                                           1294
11521 if(m .ne. 1)goto 11531                                               1294
      alm=big                                                              1294
      goto 11541                                                           1295
11531 continue                                                             1295
      alm=0.0                                                              1296
11550 do 11551 j=1,ni                                                      1296
      if(ju(j).eq.0)goto 11551                                             1296
      if(vp(j).le.0.0)goto 11551                                           1297
      alm=max(alm,abs(g(j))/vp(j))                                         1298
11551 continue                                                             1299
11552 continue                                                             1299
      alm=alf*alm/max(bta,1.0d-3)                                          1300
11541 continue                                                             1301
11501 continue                                                             1301
      dem=alm*omb                                                          1301
      ab=alm*bta                                                           1301
      rsq0=rsq                                                             1301
      jz=1                                                                 1302
11560 continue                                                             1302
11561 continue                                                             1302
      if(iz*jz.ne.0) go to 10360                                           1302
      nlp=nlp+1                                                            1302
      dlx=0.0                                                              1303
11570 do 11571 k=1,ni                                                      1303
      if(ju(k).eq.0)goto 11571                                             1304
      ak=a(k)                                                              1304
      u=g(k)+ak*xv(k)                                                      1304
      v=abs(u)-vp(k)*ab                                                    1304
      a(k)=0.0                                                             1306
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   1307 
     *em)))
      if(a(k).eq.ak)goto 11571                                             1308
      if(mm(k) .ne. 0)goto 11591                                           1308
      nin=nin+1                                                            1308
      if(nin.gt.nx)goto 11572                                              1309
11600 do 11601 j=1,ni                                                      1309
      if(ju(j).eq.0)goto 11601                                             1310
      if(mm(j) .eq. 0)goto 11621                                           1310
      c(j,nin)=c(k,mm(j))                                                  1310
      goto 11601                                                           1310
11621 continue                                                             1311
      if(j .ne. k)goto 11641                                               1311
      c(j,nin)=xv(j)                                                       1311
      goto 11601                                                           1311
11641 continue                                                             1312
      c(j,nin)=  (row_prod(j,k,ix,jx,x,w)-xm(j)*xm(k))/(xs(j)*xs(k))       1314
11601 continue                                                             1315
11602 continue                                                             1315
      mm(k)=nin                                                            1315
      ia(nin)=k                                                            1316
11591 continue                                                             1317
      del=a(k)-ak                                                          1317
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                     1318
      dlx=max(xv(k)*del**2,dlx)                                            1319
11650 do 11651 j=1,ni                                                      1319
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del                              1319
11651 continue                                                             1320
11652 continue                                                             1320
11571 continue                                                             1321
11572 continue                                                             1321
      if(dlx.lt.thr)goto 11562                                             1321
      if(nin.gt.nx)goto 11562                                              1322
      if(nlp .le. maxit)goto 11671                                         1322
      jerr=-m                                                              1322
      return                                                               1322
11671 continue                                                             1323
10360 continue                                                             1323
      iz=1                                                                 1323
      da(1:nin)=a(ia(1:nin))                                               1324
11680 continue                                                             1324
11681 continue                                                             1324
      nlp=nlp+1                                                            1324
      dlx=0.0                                                              1325
11690 do 11691 l=1,nin                                                     1325
      k=ia(l)                                                              1326
      ak=a(k)                                                              1326
      u=g(k)+ak*xv(k)                                                      1326
      v=abs(u)-vp(k)*ab                                                    1326
      a(k)=0.0                                                             1328
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   1329 
     *em)))
      if(a(k).eq.ak)goto 11691                                             1330
      del=a(k)-ak                                                          1330
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                     1331
      dlx=max(xv(k)*del**2,dlx)                                            1332
11700 do 11701 j=1,nin                                                     1332
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                                 1332
11701 continue                                                             1333
11702 continue                                                             1333
11691 continue                                                             1334
11692 continue                                                             1334
      if(dlx.lt.thr)goto 11682                                             1334
      if(nlp .le. maxit)goto 11721                                         1334
      jerr=-m                                                              1334
      return                                                               1334
11721 continue                                                             1335
      goto 11681                                                           1336
11682 continue                                                             1336
      da(1:nin)=a(ia(1:nin))-da(1:nin)                                     1337
11730 do 11731 j=1,ni                                                      1337
      if(mm(j).ne.0)goto 11731                                             1338
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))           1339
11731 continue                                                             1340
11732 continue                                                             1340
      jz=0                                                                 1341
      goto 11561                                                           1342
11562 continue                                                             1342
      if(nin .le. nx)goto 11751                                            1342
      jerr=-10000-m                                                        1342
      goto 11492                                                           1342
11751 continue                                                             1343
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                1343
      kin(m)=nin                                                           1344
      rsqo(m)=rsq                                                          1344
      almo(m)=alm                                                          1344
      lmu=m                                                                1345
      if(m.lt.mnl)goto 11491                                               1345
      if(flmin.ge.1.0)goto 11491                                           1346
      me=0                                                                 1346
11760 do 11761 j=1,nin                                                     1346
      if(ao(j,m).ne.0.0) me=me+1                                           1346
11761 continue                                                             1346
11762 continue                                                             1346
      if(me.gt.ne)goto 11492                                               1347
      if(rsq-rsq0.lt.sml*rsq)goto 11492                                    1347
      if(rsq.gt.rsqmax)goto 11492                                          1348
11491 continue                                                             1349
11492 continue                                                             1349
      deallocate(a,mm,c,da)                                                1350
      return                                                               1351
      end                                                                  1352
      subroutine spelnetn(parm,no,ni,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam,flm   1354 
     *in,ulam,  thr,isd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   1355
      double precision x(*),vp(ni),y(no),w(no),ulam(nlam),cl(2,ni)         1356
      double precision ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)            1357
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1358
      double precision, dimension (:), allocatable :: xm,xs,xv,vlam             
      integer, dimension (:), allocatable :: ju                                 
      allocate(xm(1:ni),stat=jerr)                                         1363
      if(jerr.ne.0) return                                                 1364
      allocate(xs(1:ni),stat=jerr)                                         1365
      if(jerr.ne.0) return                                                 1366
      allocate(ju(1:ni),stat=jerr)                                         1367
      if(jerr.ne.0) return                                                 1368
      allocate(xv(1:ni),stat=jerr)                                         1369
      if(jerr.ne.0) return                                                 1370
      allocate(vlam(1:nlam),stat=jerr)                                     1371
      if(jerr.ne.0) return                                                 1372
      call spchkvars(no,ni,x,ix,ju)                                        1373
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1374
      if(maxval(ju) .gt. 0)goto 11781                                      1374
      jerr=7777                                                            1374
      return                                                               1374
11781 continue                                                             1375
      call spstandard1(no,ni,x,ix,jx,y,w,ju,isd,intr,xm,xs,ym,ys,xv,jerr   1376 
     *)
      if(jerr.ne.0) return                                                 1377
      cl=cl/ys                                                             1377
      if(isd .le. 0)goto 11801                                             1377
11810 do 11811 j=1,ni                                                      1377
      cl(:,j)=cl(:,j)*xs(j)                                                1377
11811 continue                                                             1377
11812 continue                                                             1377
11801 continue                                                             1378
      if(flmin.ge.1.0) vlam=ulam/ys                                        1379
      call spelnet2(parm,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,cl,nlam,flmin,vla   1381 
     *m,thr,maxit,  xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 1382
11820 do 11821 k=1,lmu                                                     1382
      alm(k)=ys*alm(k)                                                     1382
      nk=nin(k)                                                            1383
11830 do 11831 l=1,nk                                                      1383
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                         1383
11831 continue                                                             1383
11832 continue                                                             1383
      a0(k)=0.0                                                            1384
      if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))          1385
11821 continue                                                             1386
11822 continue                                                             1386
      deallocate(xm,xs,ju,xv,vlam)                                         1387
      return                                                               1388
      end                                                                  1389
      subroutine spstandard1 (no,ni,x,ix,jx,y,w,ju,isd,intr,xm,xs,ym,ys,   1390 
     *xv,jerr)
      implicit double precision(a-h,o-z)                                   1391
      double precision x(*),y(no),w(no),xm(ni),xs(ni),xv(ni)               1392
      integer ix(*),jx(*),ju(ni)                                           1393
      w=w/sum(w)                                                           1394
      if(intr .ne. 0)goto 11851                                            1394
      ym=0.0                                                               1395
      ys=sqrt(dot_product(w,y**2)-dot_product(w,y)**2)                     1395
      y=y/ys                                                               1396
11860 do 11861 j=1,ni                                                      1396
      if(ju(j).eq.0)goto 11861                                             1396
      xm(j)=0.0                                                            1396
      jb=ix(j)                                                             1396
      je=ix(j+1)-1                                                         1397
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)                          1398
      if(isd .eq. 0)goto 11881                                             1398
      xbq=dot_product(w(jx(jb:je)),x(jb:je))**2                            1398
      vc=xv(j)-xbq                                                         1399
      xs(j)=sqrt(vc)                                                       1399
      xv(j)=1.0+xbq/vc                                                     1400
      goto 11891                                                           1401
11881 continue                                                             1401
      xs(j)=1.0                                                            1401
11891 continue                                                             1402
11871 continue                                                             1402
11861 continue                                                             1403
11862 continue                                                             1403
      return                                                               1404
11851 continue                                                             1405
11900 do 11901 j=1,ni                                                      1405
      if(ju(j).eq.0)goto 11901                                             1406
      jb=ix(j)                                                             1406
      je=ix(j+1)-1                                                         1406
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1407
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 1408
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       1409
11901 continue                                                             1410
11902 continue                                                             1410
      if(isd .ne. 0)goto 11921                                             1410
      xs=1.0                                                               1410
      goto 11931                                                           1410
11921 continue                                                             1410
      xv=1.0                                                               1410
11931 continue                                                             1411
11911 continue                                                             1411
      ym=dot_product(w,y)                                                  1411
      y=y-ym                                                               1411
      ys=sqrt(dot_product(w,y**2))                                         1411
      y=y/ys                                                               1412
      return                                                               1413
      end                                                                  1414
      subroutine spelnet2(beta,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,cl,nlam,flm   1416 
     *in,ulam,  thr,maxit,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      implicit double precision(a-h,o-z)                                   1417
      double precision y(no),w(no),x(*),vp(ni),ulam(nlam),cl(2,ni)         1418
      double precision ao(nx,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),x   1419 
     *v(ni)
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                          1420
      double precision, dimension (:), allocatable :: a,g                       
      integer, dimension (:), allocatable :: mm,iy                              
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)               1425
      allocate(a(1:ni),stat=jerr)                                          1426
      if(jerr.ne.0) return                                                 1427
      allocate(mm(1:ni),stat=jerr)                                         1428
      if(jerr.ne.0) return                                                 1429
      allocate(g(1:ni),stat=jerr)                                          1430
      if(jerr.ne.0) return                                                 1431
      allocate(iy(1:ni),stat=jerr)                                         1432
      if(jerr.ne.0) return                                                 1433
      bta=beta                                                             1433
      omb=1.0-bta                                                          1433
      alm=0.0                                                              1433
      iy=0                                                                 1435
      alf=1.0                                                              1437
      if(flmin .ge. 1.0)goto 11951                                         1437
      eqs=max(eps,flmin)                                                   1437
      alf=eqs**(1.0/(nlam-1))                                              1437
11951 continue                                                             1438
      rsq=0.0                                                              1438
      a=0.0                                                                1438
      mm=0                                                                 1438
      o=0.0                                                                1438
      nlp=0                                                                1438
      nin=nlp                                                              1438
      iz=0                                                                 1438
      mnl=min(mnlam,nlam)                                                  1439
11960 do 11961 j=1,ni                                                      1439
      if(ju(j).eq.0)goto 11961                                             1440
      jb=ix(j)                                                             1440
      je=ix(j+1)-1                                                         1441
      g(j)=abs(dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(j))    1442
11961 continue                                                             1443
11962 continue                                                             1443
11970 do 11971 m=1,nlam                                                    1443
      alm0=alm                                                             1444
      if(flmin .lt. 1.0)goto 11991                                         1444
      alm=ulam(m)                                                          1444
      goto 11981                                                           1445
11991 if(m .le. 2)goto 12001                                               1445
      alm=alm*alf                                                          1445
      goto 11981                                                           1446
12001 if(m .ne. 1)goto 12011                                               1446
      alm=big                                                              1446
      goto 12021                                                           1447
12011 continue                                                             1447
      alm0=0.0                                                             1448
12030 do 12031 j=1,ni                                                      1448
      if(ju(j).eq.0)goto 12031                                             1448
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           1448
12031 continue                                                             1449
12032 continue                                                             1449
      alm0=alm0/max(bta,1.0d-3)                                            1449
      alm=alf*alm0                                                         1450
12021 continue                                                             1451
11981 continue                                                             1451
      dem=alm*omb                                                          1451
      ab=alm*bta                                                           1451
      rsq0=rsq                                                             1451
      jz=1                                                                 1452
      tlam=bta*(2.0*alm-alm0)                                              1453
12040 do 12041 k=1,ni                                                      1453
      if(iy(k).eq.1)goto 12041                                             1453
      if(ju(k).eq.0)goto 12041                                             1454
      if(g(k).gt.tlam*vp(k)) iy(k)=1                                       1455
12041 continue                                                             1456
12042 continue                                                             1456
12050 continue                                                             1456
12051 continue                                                             1456
      if(iz*jz.ne.0) go to 10360                                           1457
11020 continue                                                             1457
      nlp=nlp+1                                                            1457
      dlx=0.0                                                              1458
12060 do 12061 k=1,ni                                                      1458
      if(iy(k).eq.0)goto 12061                                             1458
      jb=ix(k)                                                             1458
      je=ix(k+1)-1                                                         1459
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)           1460
      ak=a(k)                                                              1460
      u=gk+ak*xv(k)                                                        1460
      v=abs(u)-vp(k)*ab                                                    1460
      a(k)=0.0                                                             1462
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   1463 
     *em)))
      if(a(k).eq.ak)goto 12061                                             1464
      if(mm(k) .ne. 0)goto 12081                                           1464
      nin=nin+1                                                            1464
      if(nin.gt.nx)goto 12062                                              1465
      mm(k)=nin                                                            1465
      ia(nin)=k                                                            1466
12081 continue                                                             1467
      del=a(k)-ak                                                          1467
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1468
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)                         1469
      o=o+del*xm(k)/xs(k)                                                  1469
      dlx=max(xv(k)*del**2,dlx)                                            1470
12061 continue                                                             1471
12062 continue                                                             1471
      if(nin.gt.nx)goto 12052                                              1472
      if(dlx .ge. thr)goto 12101                                           1472
      ixx=0                                                                1473
12110 do 12111 j=1,ni                                                      1473
      if(iy(j).eq.1)goto 12111                                             1473
      if(ju(j).eq.0)goto 12111                                             1474
      jb=ix(j)                                                             1474
      je=ix(j+1)-1                                                         1475
      g(j)=abs(dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(j))    1476
      if(g(j) .le. ab*vp(j))goto 12131                                     1476
      iy(j)=1                                                              1476
      ixx=1                                                                1476
12131 continue                                                             1477
12111 continue                                                             1478
12112 continue                                                             1478
      if(ixx.eq.1) go to 11020                                             1479
      goto 12052                                                           1480
12101 continue                                                             1481
      if(nlp .le. maxit)goto 12151                                         1481
      jerr=-m                                                              1481
      return                                                               1481
12151 continue                                                             1482
10360 continue                                                             1482
      iz=1                                                                 1483
12160 continue                                                             1483
12161 continue                                                             1483
      nlp=nlp+1                                                            1483
      dlx=0.0                                                              1484
12170 do 12171 l=1,nin                                                     1484
      k=ia(l)                                                              1484
      jb=ix(k)                                                             1484
      je=ix(k+1)-1                                                         1485
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)           1486
      ak=a(k)                                                              1486
      u=gk+ak*xv(k)                                                        1486
      v=abs(u)-vp(k)*ab                                                    1486
      a(k)=0.0                                                             1488
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   1489 
     *em)))
      if(a(k).eq.ak)goto 12171                                             1490
      del=a(k)-ak                                                          1490
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1491
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)                         1492
      o=o+del*xm(k)/xs(k)                                                  1492
      dlx=max(xv(k)*del**2,dlx)                                            1493
12171 continue                                                             1494
12172 continue                                                             1494
      if(dlx.lt.thr)goto 12162                                             1494
      if(nlp .le. maxit)goto 12191                                         1494
      jerr=-m                                                              1494
      return                                                               1494
12191 continue                                                             1495
      goto 12161                                                           1496
12162 continue                                                             1496
      jz=0                                                                 1497
      goto 12051                                                           1498
12052 continue                                                             1498
      if(nin .le. nx)goto 12211                                            1498
      jerr=-10000-m                                                        1498
      goto 11972                                                           1498
12211 continue                                                             1499
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                1499
      kin(m)=nin                                                           1500
      rsqo(m)=rsq                                                          1500
      almo(m)=alm                                                          1500
      lmu=m                                                                1501
      if(m.lt.mnl)goto 11971                                               1501
      if(flmin.ge.1.0)goto 11971                                           1502
      me=0                                                                 1502
12220 do 12221 j=1,nin                                                     1502
      if(ao(j,m).ne.0.0) me=me+1                                           1502
12221 continue                                                             1502
12222 continue                                                             1502
      if(me.gt.ne)goto 11972                                               1503
      if(rsq-rsq0.lt.sml*rsq)goto 11972                                    1503
      if(rsq.gt.rsqmax)goto 11972                                          1504
11971 continue                                                             1505
11972 continue                                                             1505
      deallocate(a,mm,g,iy)                                                1506
      return                                                               1507
      end                                                                  1508
      subroutine spchkvars(no,ni,x,ix,ju)                                  1509
      implicit double precision(a-h,o-z)                                   1510
      double precision x(*)                                                1510
      integer ix(*),ju(ni)                                                 1511
12230 do 12231 j=1,ni                                                      1511
      ju(j)=0                                                              1511
      jb=ix(j)                                                             1511
      nj=ix(j+1)-jb                                                        1511
      if(nj.eq.0)goto 12231                                                1512
      je=ix(j+1)-1                                                         1513
      if(nj .ge. no)goto 12251                                             1513
12260 do 12261 i=jb,je                                                     1513
      if(x(i).eq.0.0)goto 12261                                            1513
      ju(j)=1                                                              1513
      goto 12262                                                           1513
12261 continue                                                             1513
12262 continue                                                             1513
      goto 12271                                                           1514
12251 continue                                                             1514
      t=x(jb)                                                              1514
12280 do 12281 i=jb+1,je                                                   1514
      if(x(i).eq.t)goto 12281                                              1514
      ju(j)=1                                                              1514
      goto 12282                                                           1514
12281 continue                                                             1514
12282 continue                                                             1514
12271 continue                                                             1515
12241 continue                                                             1515
12231 continue                                                             1516
12232 continue                                                             1516
      return                                                               1517
      end                                                                  1518
      subroutine cmodval(a0,ca,ia,nin,x,ix,jx,n,f)                         1519
      implicit double precision(a-h,o-z)                                   1520
      double precision ca(*),x(*),f(n)                                     1520
      integer ia(*),ix(*),jx(*)                                            1521
      f=a0                                                                 1522
12290 do 12291 j=1,nin                                                     1522
      k=ia(j)                                                              1522
      kb=ix(k)                                                             1522
      ke=ix(k+1)-1                                                         1523
      f(jx(kb:ke))=f(jx(kb:ke))+ca(j)*x(kb:ke)                             1524
12291 continue                                                             1525
12292 continue                                                             1525
      return                                                               1526
      end                                                                  1527
      function row_prod(i,j,ia,ja,ra,w)                                    1528
      implicit double precision(a-h,o-z)                                   1529
      integer ia(*),ja(*)                                                  1529
      double precision ra(*),w(*)                                          1530
      row_prod=dot(ra(ia(i)),ra(ia(j)),ja(ia(i)),ja(ia(j)),  ia(i+1)-ia(   1532 
     *i),ia(j+1)-ia(j),w)
      return                                                               1533
      end                                                                  1534
      function dot(x,y,mx,my,nx,ny,w)                                      1535
      implicit double precision(a-h,o-z)                                   1536
      double precision x(*),y(*),w(*)                                      1536
      integer mx(*),my(*)                                                  1537
      i=1                                                                  1537
      j=i                                                                  1537
      s=0.0                                                                1538
12300 continue                                                             1538
12301 continue                                                             1538
12310 continue                                                             1539
12311 if(mx(i).ge.my(j))goto 12312                                         1539
      i=i+1                                                                1539
      if(i.gt.nx) go to 12320                                              1539
      goto 12311                                                           1540
12312 continue                                                             1540
      if(mx(i).eq.my(j)) go to 12330                                       1541
12340 continue                                                             1541
12341 if(my(j).ge.mx(i))goto 12342                                         1541
      j=j+1                                                                1541
      if(j.gt.ny) go to 12320                                              1541
      goto 12341                                                           1542
12342 continue                                                             1542
      if(mx(i).eq.my(j)) go to 12330                                       1542
      goto 12301                                                           1543
12330 continue                                                             1543
      s=s+w(mx(i))*x(i)*y(j)                                               1544
      i=i+1                                                                1544
      if(i.gt.nx)goto 12302                                                1544
      j=j+1                                                                1544
      if(j.gt.ny)goto 12302                                                1545
      goto 12301                                                           1546
12302 continue                                                             1546
12320 continue                                                             1546
      dot=s                                                                1547
      return                                                               1548
      end                                                                  1549
      subroutine lognet (parm,no,ni,nc,x,y,g,jd,vp,cl,ne,nx,nlam,flmin,u   1551 
     *lam,thr,  isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,je
     *rr)
      implicit double precision(a-h,o-z)                                   1552
      double precision x(no,ni),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nla   1553 
     *m)
      double precision ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl   1554 
     *(2,ni)
      integer jd(*),ia(nx),nin(nlam)                                       1555
      double precision, dimension (:), allocatable :: xm,xs,ww,vq,xv            
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 12361                                    1559
      jerr=10000                                                           1559
      return                                                               1559
12361 continue                                                             1560
      allocate(ww(1:no),stat=jerr)                                         1561
      if(jerr.ne.0) return                                                 1562
      allocate(ju(1:ni),stat=jerr)                                         1563
      if(jerr.ne.0) return                                                 1564
      allocate(vq(1:ni),stat=jerr)                                         1565
      if(jerr.ne.0) return                                                 1566
      allocate(xm(1:ni),stat=jerr)                                         1567
      if(jerr.ne.0) return                                                 1568
      if(kopt .ne. 2)goto 12381                                            1568
      allocate(xv(1:ni),stat=jerr)                                         1568
      if(jerr.ne.0) return                                                 1568
12381 continue                                                             1569
      if(isd .le. 0)goto 12401                                             1569
      allocate(xs(1:ni),stat=jerr)                                         1569
      if(jerr.ne.0) return                                                 1569
12401 continue                                                             1571
      call chkvars(no,ni,x,ju)                                             1572
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1573
      if(maxval(ju) .gt. 0)goto 12421                                      1573
      jerr=7777                                                            1573
      return                                                               1573
12421 continue                                                             1574
      vq=max(0d0,vp)                                                       1574
      vq=vq*ni/sum(vq)                                                     1575
12430 do 12431 i=1,no                                                      1575
      ww(i)=sum(y(i,:))                                                    1575
      if(ww(i).gt.0.0) y(i,:)=y(i,:)/ww(i)                                 1575
12431 continue                                                             1576
12432 continue                                                             1576
      sw=sum(ww)                                                           1576
      ww=ww/sw                                                             1577
      if(nc .ne. 1)goto 12451                                              1577
      call lstandard1(no,ni,x,ww,ju,isd,intr,xm,xs)                        1578
      if(isd .le. 0)goto 12471                                             1578
12480 do 12481 j=1,ni                                                      1578
      cl(:,j)=cl(:,j)*xs(j)                                                1578
12481 continue                                                             1578
12482 continue                                                             1578
12471 continue                                                             1579
      call lognet2n(parm,no,ni,x,y(:,1),g(:,1),ww,ju,vq,cl,ne,nx,nlam,fl   1581 
     *min,ulam,  thr,isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,n
     *lp,jerr)
      goto 12441                                                           1582
12451 if(kopt .ne. 2)goto 12491                                            1582
      call multlstandard1(no,ni,x,ww,ju,isd,intr,xm,xs,xv)                 1583
      if(isd .le. 0)goto 12511                                             1583
12520 do 12521 j=1,ni                                                      1583
      cl(:,j)=cl(:,j)*xs(j)                                                1583
12521 continue                                                             1583
12522 continue                                                             1583
12511 continue                                                             1584
      call multlognetn(parm,no,ni,nc,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,   1586 
     *ulam,thr,  intr,maxit,xv,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      goto 12531                                                           1587
12491 continue                                                             1587
      call lstandard1(no,ni,x,ww,ju,isd,intr,xm,xs)                        1588
      if(isd .le. 0)goto 12551                                             1588
12560 do 12561 j=1,ni                                                      1588
      cl(:,j)=cl(:,j)*xs(j)                                                1588
12561 continue                                                             1588
12562 continue                                                             1588
12551 continue                                                             1589
      call lognetn(parm,no,ni,nc,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam   1591 
     *,thr,  isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
12531 continue                                                             1592
12441 continue                                                             1592
      if(jerr.gt.0) return                                                 1592
      dev0=2.0*sw*dev0                                                     1593
12570 do 12571 k=1,lmu                                                     1593
      nk=nin(k)                                                            1594
12580 do 12581 ic=1,nc                                                     1594
      if(isd .le. 0)goto 12601                                             1594
12610 do 12611 l=1,nk                                                      1594
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                      1594
12611 continue                                                             1594
12612 continue                                                             1594
12601 continue                                                             1595
      if(intr .ne. 0)goto 12631                                            1595
      a0(ic,k)=0.0                                                         1595
      goto 12641                                                           1596
12631 continue                                                             1596
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))            1596
12641 continue                                                             1597
12621 continue                                                             1597
12581 continue                                                             1598
12582 continue                                                             1598
12571 continue                                                             1599
12572 continue                                                             1599
      deallocate(ww,ju,vq,xm)                                              1599
      if(isd.gt.0) deallocate(xs)                                          1600
      if(kopt.eq.2) deallocate(xv)                                         1601
      return                                                               1602
      end                                                                  1603
      subroutine lstandard1 (no,ni,x,w,ju,isd,intr,xm,xs)                  1604
      implicit double precision(a-h,o-z)                                   1605
      double precision x(no,ni),w(no),xm(ni),xs(ni)                        1605
      integer ju(ni)                                                       1606
      if(intr .ne. 0)goto 12661                                            1607
12670 do 12671 j=1,ni                                                      1607
      if(ju(j).eq.0)goto 12671                                             1607
      xm(j)=0.0                                                            1608
      if(isd .eq. 0)goto 12691                                             1608
      vc=dot_product(w,x(:,j)**2)-dot_product(w,x(:,j))**2                 1609
      xs(j)=sqrt(vc)                                                       1609
      x(:,j)=x(:,j)/xs(j)                                                  1610
12691 continue                                                             1611
12671 continue                                                             1612
12672 continue                                                             1612
      return                                                               1613
12661 continue                                                             1614
12700 do 12701 j=1,ni                                                      1614
      if(ju(j).eq.0)goto 12701                                             1615
      xm(j)=dot_product(w,x(:,j))                                          1615
      x(:,j)=x(:,j)-xm(j)                                                  1616
      if(isd .le. 0)goto 12721                                             1616
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                                 1616
      x(:,j)=x(:,j)/xs(j)                                                  1616
12721 continue                                                             1617
12701 continue                                                             1618
12702 continue                                                             1618
      return                                                               1619
      end                                                                  1620
      subroutine multlstandard1 (no,ni,x,w,ju,isd,intr,xm,xs,xv)           1621
      implicit double precision(a-h,o-z)                                   1622
      double precision x(no,ni),w(no),xm(ni),xs(ni),xv(ni)                 1622
      integer ju(ni)                                                       1623
      if(intr .ne. 0)goto 12741                                            1624
12750 do 12751 j=1,ni                                                      1624
      if(ju(j).eq.0)goto 12751                                             1624
      xm(j)=0.0                                                            1625
      xv(j)=dot_product(w,x(:,j)**2)                                       1626
      if(isd .eq. 0)goto 12771                                             1626
      xbq=dot_product(w,x(:,j))**2                                         1626
      vc=xv(j)-xbq                                                         1627
      xs(j)=sqrt(vc)                                                       1627
      x(:,j)=x(:,j)/xs(j)                                                  1627
      xv(j)=1.0+xbq/vc                                                     1628
12771 continue                                                             1629
12751 continue                                                             1630
12752 continue                                                             1630
      return                                                               1631
12741 continue                                                             1632
12780 do 12781 j=1,ni                                                      1632
      if(ju(j).eq.0)goto 12781                                             1633
      xm(j)=dot_product(w,x(:,j))                                          1633
      x(:,j)=x(:,j)-xm(j)                                                  1634
      xv(j)=dot_product(w,x(:,j)**2)                                       1635
      if(isd .le. 0)goto 12801                                             1635
      xs(j)=sqrt(xv(j))                                                    1635
      x(:,j)=x(:,j)/xs(j)                                                  1635
      xv(j)=1.0                                                            1635
12801 continue                                                             1636
12781 continue                                                             1637
12782 continue                                                             1637
      return                                                               1638
      end                                                                  1639
      subroutine lognet2n(parm,no,ni,x,y,g,w,ju,vp,cl,ne,nx,nlam,flmin,u   1641 
     *lam,shri,  isd,intr,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jer
     *r)
      implicit double precision(a-h,o-z)                                   1642
      double precision x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam),cl(2   1643 
     *,ni)
      double precision a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)             1644
      integer ju(ni),m(nx),kin(nlam)                                       1645
      double precision, dimension (:), allocatable :: b,bs,v,r,xv,q,ga          
      integer, dimension (:), allocatable :: mm,ixx                             
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               1650
      allocate(b(0:ni),stat=jerr)                                          1651
      if(jerr.ne.0) return                                                 1652
      allocate(xv(1:ni),stat=jerr)                                         1653
      if(jerr.ne.0) return                                                 1654
      allocate(ga(1:ni),stat=jerr)                                         1655
      if(jerr.ne.0) return                                                 1656
      allocate(bs(0:ni),stat=jerr)                                         1657
      if(jerr.ne.0) return                                                 1658
      allocate(mm(1:ni),stat=jerr)                                         1659
      if(jerr.ne.0) return                                                 1660
      allocate(ixx(1:ni),stat=jerr)                                        1661
      if(jerr.ne.0) return                                                 1662
      allocate(r(1:no),stat=jerr)                                          1663
      if(jerr.ne.0) return                                                 1664
      allocate(v(1:no),stat=jerr)                                          1665
      if(jerr.ne.0) return                                                 1666
      allocate(q(1:no),stat=jerr)                                          1667
      if(jerr.ne.0) return                                                 1668
      fmax=log(1.0/pmin-1.0)                                               1668
      fmin=-fmax                                                           1668
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                      1669
      bta=parm                                                             1669
      omb=1.0-bta                                                          1670
      q0=dot_product(w,y)                                                  1670
      if(q0 .gt. pmin)goto 12821                                           1670
      jerr=8001                                                            1670
      return                                                               1670
12821 continue                                                             1671
      if(q0 .lt. 1.0-pmin)goto 12841                                       1671
      jerr=9001                                                            1671
      return                                                               1671
12841 continue                                                             1672
      if(intr.eq.0.0) q0=0.5                                               1673
      ixx=0                                                                1673
      al=0.0                                                               1673
      bz=0.0                                                               1673
      if(intr.ne.0) bz=log(q0/(1.0-q0))                                    1674
      if(nonzero(no,g) .ne. 0)goto 12861                                   1674
      vi=q0*(1.0-q0)                                                       1674
      b(0)=bz                                                              1674
      v=vi*w                                                               1675
      r=w*(y-q0)                                                           1675
      q=q0                                                                 1675
      xmz=vi                                                               1675
      dev1=-(bz*q0+log(1.0-q0))                                            1676
      goto 12871                                                           1677
12861 continue                                                             1677
      b(0)=0.0                                                             1678
      if(intr .eq. 0)goto 12891                                            1678
      b(0)=azero(no,y,g,w,jerr)                                            1678
      if(jerr.ne.0) return                                                 1678
12891 continue                                                             1679
      q=1.0/(1.0+exp(-b(0)-g))                                             1679
      v=w*q*(1.0-q)                                                        1679
      r=w*(y-q)                                                            1679
      xmz=sum(v)                                                           1680
      dev1=-(b(0)*q0+dot_product(w,y*g+log(1.0-q)))                        1681
12871 continue                                                             1682
12851 continue                                                             1682
      if(kopt .le. 0)goto 12911                                            1683
      if(isd .le. 0 .or. intr .eq. 0)goto 12931                            1683
      xv=0.25                                                              1683
      goto 12941                                                           1684
12931 continue                                                             1684
12950 do 12951 j=1,ni                                                      1684
      if(ju(j).ne.0) xv(j)=0.25*dot_product(w,x(:,j)**2)                   1684
12951 continue                                                             1684
12952 continue                                                             1684
12941 continue                                                             1685
12921 continue                                                             1685
12911 continue                                                             1686
      dev0=dev1                                                            1687
12960 do 12961 i=1,no                                                      1687
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                        1688
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))              1689
12961 continue                                                             1691
12962 continue                                                             1691
      alf=1.0                                                              1693
      if(flmin .ge. 1.0)goto 12981                                         1693
      eqs=max(eps,flmin)                                                   1693
      alf=eqs**(1.0/(nlam-1))                                              1693
12981 continue                                                             1694
      m=0                                                                  1694
      mm=0                                                                 1694
      nlp=0                                                                1694
      nin=nlp                                                              1694
      mnl=min(mnlam,nlam)                                                  1694
      bs=0.0                                                               1694
      b(1:ni)=0.0                                                          1695
      shr=shri*dev0                                                        1696
12990 do 12991 j=1,ni                                                      1696
      if(ju(j).eq.0)goto 12991                                             1696
      ga(j)=abs(dot_product(r,x(:,j)))                                     1696
12991 continue                                                             1697
12992 continue                                                             1697
13000 do 13001 ilm=1,nlam                                                  1697
      al0=al                                                               1698
      if(flmin .lt. 1.0)goto 13021                                         1698
      al=ulam(ilm)                                                         1698
      goto 13011                                                           1699
13021 if(ilm .le. 2)goto 13031                                             1699
      al=al*alf                                                            1699
      goto 13011                                                           1700
13031 if(ilm .ne. 1)goto 13041                                             1700
      al=big                                                               1700
      goto 13051                                                           1701
13041 continue                                                             1701
      al0=0.0                                                              1702
13060 do 13061 j=1,ni                                                      1702
      if(ju(j).eq.0)goto 13061                                             1702
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            1702
13061 continue                                                             1703
13062 continue                                                             1703
      al0=al0/max(bta,1.0d-3)                                              1703
      al=alf*al0                                                           1704
13051 continue                                                             1705
13011 continue                                                             1705
      al2=al*omb                                                           1705
      al1=al*bta                                                           1705
      tlam=bta*(2.0*al-al0)                                                1706
13070 do 13071 k=1,ni                                                      1706
      if(ixx(k).eq.1)goto 13071                                            1706
      if(ju(k).eq.0)goto 13071                                             1707
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     1708
13071 continue                                                             1709
13072 continue                                                             1709
11020 continue                                                             1710
13080 continue                                                             1710
13081 continue                                                             1710
      bs(0)=b(0)                                                           1710
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                1711
      if(kopt .ne. 0)goto 13101                                            1712
13110 do 13111 j=1,ni                                                      1712
      if(ixx(j).gt.0) xv(j)=dot_product(v,x(:,j)**2)                       1712
13111 continue                                                             1713
13112 continue                                                             1713
13101 continue                                                             1714
13120 continue                                                             1714
13121 continue                                                             1714
      nlp=nlp+1                                                            1714
      dlx=0.0                                                              1715
13130 do 13131 k=1,ni                                                      1715
      if(ixx(k).eq.0)goto 13131                                            1716
      bk=b(k)                                                              1716
      gk=dot_product(r,x(:,k))                                             1717
      u=gk+xv(k)*b(k)                                                      1717
      au=abs(u)-vp(k)*al1                                                  1718
      if(au .gt. 0.0)goto 13151                                            1718
      b(k)=0.0                                                             1718
      goto 13161                                                           1719
13151 continue                                                             1720
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))          1721
13161 continue                                                             1722
13141 continue                                                             1722
      d=b(k)-bk                                                            1722
      if(abs(d).le.0.0)goto 13131                                          1722
      dlx=max(dlx,xv(k)*d**2)                                              1723
      r=r-d*v*x(:,k)                                                       1724
      if(mm(k) .ne. 0)goto 13181                                           1724
      nin=nin+1                                                            1724
      if(nin.gt.nx)goto 13132                                              1725
      mm(k)=nin                                                            1725
      m(nin)=k                                                             1726
13181 continue                                                             1727
13131 continue                                                             1728
13132 continue                                                             1728
      if(nin.gt.nx)goto 13122                                              1729
      d=0.0                                                                1729
      if(intr.ne.0) d=sum(r)/xmz                                           1730
      if(d .eq. 0.0)goto 13201                                             1730
      b(0)=b(0)+d                                                          1730
      dlx=max(dlx,xmz*d**2)                                                1730
      r=r-d*v                                                              1730
13201 continue                                                             1731
      if(dlx.lt.shr)goto 13122                                             1731
      if(nlp .le. maxit)goto 13221                                         1731
      jerr=-ilm                                                            1731
      return                                                               1731
13221 continue                                                             1732
13230 continue                                                             1732
13231 continue                                                             1732
      nlp=nlp+1                                                            1732
      dlx=0.0                                                              1733
13240 do 13241 l=1,nin                                                     1733
      k=m(l)                                                               1733
      bk=b(k)                                                              1734
      gk=dot_product(r,x(:,k))                                             1735
      u=gk+xv(k)*b(k)                                                      1735
      au=abs(u)-vp(k)*al1                                                  1736
      if(au .gt. 0.0)goto 13261                                            1736
      b(k)=0.0                                                             1736
      goto 13271                                                           1737
13261 continue                                                             1738
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))          1739
13271 continue                                                             1740
13251 continue                                                             1740
      d=b(k)-bk                                                            1740
      if(abs(d).le.0.0)goto 13241                                          1740
      dlx=max(dlx,xv(k)*d**2)                                              1741
      r=r-d*v*x(:,k)                                                       1742
13241 continue                                                             1743
13242 continue                                                             1743
      d=0.0                                                                1743
      if(intr.ne.0) d=sum(r)/xmz                                           1744
      if(d .eq. 0.0)goto 13291                                             1744
      b(0)=b(0)+d                                                          1744
      dlx=max(dlx,xmz*d**2)                                                1744
      r=r-d*v                                                              1744
13291 continue                                                             1745
      if(dlx.lt.shr)goto 13232                                             1745
      if(nlp .le. maxit)goto 13311                                         1745
      jerr=-ilm                                                            1745
      return                                                               1745
13311 continue                                                             1746
      goto 13231                                                           1747
13232 continue                                                             1747
      goto 13121                                                           1748
13122 continue                                                             1748
      if(nin.gt.nx)goto 13082                                              1749
13320 do 13321 i=1,no                                                      1749
      fi=b(0)+g(i)                                                         1750
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin)),x(i,m(1:nin)))            1751
      if(fi .ge. fmin)goto 13341                                           1751
      q(i)=0.0                                                             1751
      goto 13331                                                           1751
13341 if(fi .le. fmax)goto 13351                                           1751
      q(i)=1.0                                                             1751
      goto 13361                                                           1752
13351 continue                                                             1752
      q(i)=1.0/(1.0+exp(-fi))                                              1752
13361 continue                                                             1753
13331 continue                                                             1753
13321 continue                                                             1754
13322 continue                                                             1754
      v=w*q*(1.0-q)                                                        1754
      xmz=sum(v)                                                           1754
      if(xmz.le.vmin)goto 13082                                            1754
      r=w*(y-q)                                                            1755
      if(xmz*(b(0)-bs(0))**2 .ge. shr)goto 13381                           1755
      ix=0                                                                 1756
13390 do 13391 j=1,nin                                                     1756
      k=m(j)                                                               1757
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 13391                           1757
      ix=1                                                                 1757
      goto 13392                                                           1758
13391 continue                                                             1759
13392 continue                                                             1759
      if(ix .ne. 0)goto 13411                                              1760
13420 do 13421 k=1,ni                                                      1760
      if(ixx(k).eq.1)goto 13421                                            1760
      if(ju(k).eq.0)goto 13421                                             1761
      ga(k)=abs(dot_product(r,x(:,k)))                                     1762
      if(ga(k) .le. al1*vp(k))goto 13441                                   1762
      ixx(k)=1                                                             1762
      ix=1                                                                 1762
13441 continue                                                             1763
13421 continue                                                             1764
13422 continue                                                             1764
      if(ix.eq.1) go to 11020                                              1765
      goto 13082                                                           1766
13411 continue                                                             1767
13381 continue                                                             1768
      goto 13081                                                           1769
13082 continue                                                             1769
      if(nin .le. nx)goto 13461                                            1769
      jerr=-10000-ilm                                                      1769
      goto 13002                                                           1769
13461 continue                                                             1770
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                1770
      kin(ilm)=nin                                                         1771
      a0(ilm)=b(0)                                                         1771
      alm(ilm)=al                                                          1771
      lmu=ilm                                                              1772
      devi=dev2(no,w,y,q,pmin)                                             1773
      dev(ilm)=(dev1-devi)/dev0                                            1773
      if(xmz.le.vmin)goto 13002                                            1774
      if(ilm.lt.mnl)goto 13001                                             1774
      if(flmin.ge.1.0)goto 13001                                           1775
      me=0                                                                 1775
13470 do 13471 j=1,nin                                                     1775
      if(a(j,ilm).ne.0.0) me=me+1                                          1775
13471 continue                                                             1775
13472 continue                                                             1775
      if(me.gt.ne)goto 13002                                               1776
      if(dev(ilm).gt.devmax)goto 13002                                     1776
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 13002                             1777
13001 continue                                                             1778
13002 continue                                                             1778
      g=log(q/(1.0-q))                                                     1779
      deallocate(b,bs,v,r,xv,q,mm,ga,ixx)                                  1780
      return                                                               1781
      end                                                                  1782
      function dev2(n,w,y,p,pmin)                                          1783
      implicit double precision(a-h,o-z)                                   1784
      double precision w(n),y(n),p(n)                                      1785
      pmax=1.0-pmin                                                        1785
      s=0.0                                                                1786
13480 do 13481 i=1,n                                                       1786
      pi=min(max(pmin,p(i)),pmax)                                          1787
      s=s-w(i)*(y(i)*log(pi)+(1.0-y(i))*log(1.0-pi))                       1788
13481 continue                                                             1789
13482 continue                                                             1789
      dev2=s                                                               1790
      return                                                               1791
      end                                                                  1792
      function azero(n,y,g,q,jerr)                                         1793
      implicit double precision(a-h,o-z)                                   1794
      parameter(eps=1.0d-7)                                                1795
      double precision y(n),g(n),q(n)                                      1796
      double precision, dimension (:), allocatable :: e,p,w                     
      azero = 0.0                                                          1800
      allocate(e(1:n),stat=jerr)                                           1801
      if(jerr.ne.0) return                                                 1802
      allocate(p(1:n),stat=jerr)                                           1803
      if(jerr.ne.0) return                                                 1804
      allocate(w(1:n),stat=jerr)                                           1805
      if(jerr.ne.0) return                                                 1806
      az=0.0                                                               1806
      e=exp(-g)                                                            1806
      qy=dot_product(q,y)                                                  1806
      p=1.0/(1.0+e)                                                        1807
13490 continue                                                             1807
13491 continue                                                             1807
      w=q*p*(1.0-p)                                                        1808
      d=(qy-dot_product(q,p))/sum(w)                                       1808
      az=az+d                                                              1808
      if(abs(d).lt.eps)goto 13492                                          1809
      ea0=exp(-az)                                                         1809
      p=1.0/(1.0+ea0*e)                                                    1810
      goto 13491                                                           1811
13492 continue                                                             1811
      azero=az                                                             1812
      deallocate(e,p,w)                                                    1813
      return                                                               1814
      end                                                                  1815
      subroutine lognetn(parm,no,ni,nc,x,y,g,w,ju,vp,cl,ne,nx,nlam,flmin   1817 
     *,ulam,shri,  isd,intr,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,j
     *err)
      implicit double precision(a-h,o-z)                                   1818
      double precision x(no,ni),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam   1819 
     *)
      double precision a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(   1820 
     *2,ni)
      integer ju(ni),m(nx),kin(nlam)                                       1821
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
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               1836
      exmn=-exmx                                                           1837
      allocate(r(1:no),stat=jerr)                                          1838
      if(jerr.ne.0) return                                                 1839
      allocate(v(1:no),stat=jerr)                                          1840
      if(jerr.ne.0) return                                                 1841
      allocate(mm(1:ni),stat=jerr)                                         1842
      if(jerr.ne.0) return                                                 1843
      allocate(is(1:max(nc,ni)),stat=jerr)                                 1844
      if(jerr.ne.0) return                                                 1845
      allocate(sxp(1:no),stat=jerr)                                        1846
      if(jerr.ne.0) return                                                 1847
      allocate(sxpl(1:no),stat=jerr)                                       1848
      if(jerr.ne.0) return                                                 1849
      allocate(di(1:no),stat=jerr)                                         1850
      if(jerr.ne.0) return                                                 1851
      allocate(ga(1:ni),stat=jerr)                                         1852
      if(jerr.ne.0) return                                                 1853
      allocate(ixx(1:ni),stat=jerr)                                        1854
      if(jerr.ne.0) return                                                 1855
      pmax=1.0-pmin                                                        1855
      emin=pmin/pmax                                                       1855
      emax=1.0/emin                                                        1856
      pfm=(1.0+pmin)*pmin                                                  1856
      pfx=(1.0-pmin)*pmax                                                  1856
      vmin=pfm*pmax                                                        1857
      bta=parm                                                             1857
      omb=1.0-bta                                                          1857
      dev1=0.0                                                             1857
      dev0=0.0                                                             1858
13500 do 13501 ic=1,nc                                                     1858
      q0=dot_product(w,y(:,ic))                                            1859
      if(q0 .gt. pmin)goto 13521                                           1859
      jerr =8000+ic                                                        1859
      return                                                               1859
13521 continue                                                             1860
      if(q0 .lt. 1.0-pmin)goto 13541                                       1860
      jerr =9000+ic                                                        1860
      return                                                               1860
13541 continue                                                             1861
      if(intr .ne. 0)goto 13561                                            1861
      q0=1.0/nc                                                            1861
      b(0,ic)=0.0                                                          1861
      goto 13571                                                           1862
13561 continue                                                             1862
      b(0,ic)=log(q0)                                                      1862
      dev1=dev1-q0*b(0,ic)                                                 1862
13571 continue                                                             1863
13551 continue                                                             1863
      b(1:ni,ic)=0.0                                                       1864
13501 continue                                                             1865
13502 continue                                                             1865
      if(intr.eq.0) dev1=log(float(nc))                                    1865
      ixx=0                                                                1865
      al=0.0                                                               1866
      if(nonzero(no*nc,g) .ne. 0)goto 13591                                1867
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         1867
      sxp=0.0                                                              1868
13600 do 13601 ic=1,nc                                                     1868
      q(:,ic)=exp(b(0,ic))                                                 1868
      sxp=sxp+q(:,ic)                                                      1868
13601 continue                                                             1869
13602 continue                                                             1869
      goto 13611                                                           1870
13591 continue                                                             1870
13620 do 13621 i=1,no                                                      1870
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1870
13621 continue                                                             1870
13622 continue                                                             1870
      sxp=0.0                                                              1871
      if(intr .ne. 0)goto 13641                                            1871
      b(0,:)=0.0                                                           1871
      goto 13651                                                           1872
13641 continue                                                             1872
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 1872
      if(jerr.ne.0) return                                                 1872
13651 continue                                                             1873
13631 continue                                                             1873
      dev1=0.0                                                             1874
13660 do 13661 ic=1,nc                                                     1874
      q(:,ic)=b(0,ic)+g(:,ic)                                              1875
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             1876
      q(:,ic)=exp(q(:,ic))                                                 1876
      sxp=sxp+q(:,ic)                                                      1877
13661 continue                                                             1878
13662 continue                                                             1878
      sxpl=w*log(sxp)                                                      1878
13670 do 13671 ic=1,nc                                                     1878
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  1878
13671 continue                                                             1879
13672 continue                                                             1879
13611 continue                                                             1880
13581 continue                                                             1880
13680 do 13681 ic=1,nc                                                     1880
13690 do 13691 i=1,no                                                      1880
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               1880
13691 continue                                                             1880
13692 continue                                                             1880
13681 continue                                                             1881
13682 continue                                                             1881
      dev0=dev0+dev1                                                       1882
      if(kopt .le. 0)goto 13711                                            1883
      if(isd .le. 0 .or. intr .eq. 0)goto 13731                            1883
      xv=0.25                                                              1883
      goto 13741                                                           1884
13731 continue                                                             1884
13750 do 13751 j=1,ni                                                      1884
      if(ju(j).ne.0) xv(j,:)=0.25*dot_product(w,x(:,j)**2)                 1884
13751 continue                                                             1884
13752 continue                                                             1884
13741 continue                                                             1885
13721 continue                                                             1885
13711 continue                                                             1887
      alf=1.0                                                              1889
      if(flmin .ge. 1.0)goto 13771                                         1889
      eqs=max(eps,flmin)                                                   1889
      alf=eqs**(1.0/(nlam-1))                                              1889
13771 continue                                                             1890
      m=0                                                                  1890
      mm=0                                                                 1890
      nin=0                                                                1890
      nlp=0                                                                1890
      mnl=min(mnlam,nlam)                                                  1890
      bs=0.0                                                               1890
      shr=shri*dev0                                                        1891
      ga=0.0                                                               1892
13780 do 13781 ic=1,nc                                                     1892
      r=w*(y(:,ic)-q(:,ic)/sxp)                                            1893
13790 do 13791 j=1,ni                                                      1893
      if(ju(j).ne.0) ga(j)=max(ga(j),abs(dot_product(r,x(:,j))))           1893
13791 continue                                                             1894
13792 continue                                                             1894
13781 continue                                                             1895
13782 continue                                                             1895
13800 do 13801 ilm=1,nlam                                                  1895
      al0=al                                                               1896
      if(flmin .lt. 1.0)goto 13821                                         1896
      al=ulam(ilm)                                                         1896
      goto 13811                                                           1897
13821 if(ilm .le. 2)goto 13831                                             1897
      al=al*alf                                                            1897
      goto 13811                                                           1898
13831 if(ilm .ne. 1)goto 13841                                             1898
      al=big                                                               1898
      goto 13851                                                           1899
13841 continue                                                             1899
      al0=0.0                                                              1900
13860 do 13861 j=1,ni                                                      1900
      if(ju(j).eq.0)goto 13861                                             1900
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            1900
13861 continue                                                             1901
13862 continue                                                             1901
      al0=al0/max(bta,1.0d-3)                                              1901
      al=alf*al0                                                           1902
13851 continue                                                             1903
13811 continue                                                             1903
      al2=al*omb                                                           1903
      al1=al*bta                                                           1903
      tlam=bta*(2.0*al-al0)                                                1904
13870 do 13871 k=1,ni                                                      1904
      if(ixx(k).eq.1)goto 13871                                            1904
      if(ju(k).eq.0)goto 13871                                             1905
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     1906
13871 continue                                                             1907
13872 continue                                                             1907
11020 continue                                                             1908
13880 continue                                                             1908
13881 continue                                                             1908
      ix=0                                                                 1908
      jx=ix                                                                1908
      ig=0                                                                 1909
13890 do 13891 ic=1,nc                                                     1909
      bs(0,ic)=b(0,ic)                                                     1910
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          1911
      xmz=0.0                                                              1912
13900 do 13901 i=1,no                                                      1912
      pic=q(i,ic)/sxp(i)                                                   1913
      if(pic .ge. pfm)goto 13921                                           1913
      pic=0.0                                                              1913
      v(i)=0.0                                                             1913
      goto 13911                                                           1914
13921 if(pic .le. pfx)goto 13931                                           1914
      pic=1.0                                                              1914
      v(i)=0.0                                                             1914
      goto 13941                                                           1915
13931 continue                                                             1915
      v(i)=w(i)*pic*(1.0-pic)                                              1915
      xmz=xmz+v(i)                                                         1915
13941 continue                                                             1916
13911 continue                                                             1916
      r(i)=w(i)*(y(i,ic)-pic)                                              1917
13901 continue                                                             1918
13902 continue                                                             1918
      if(xmz.le.vmin)goto 13891                                            1918
      ig=1                                                                 1919
      if(kopt .ne. 0)goto 13961                                            1920
13970 do 13971 j=1,ni                                                      1920
      if(ixx(j).gt.0) xv(j,ic)=dot_product(v,x(:,j)**2)                    1920
13971 continue                                                             1921
13972 continue                                                             1921
13961 continue                                                             1922
13980 continue                                                             1922
13981 continue                                                             1922
      nlp=nlp+1                                                            1922
      dlx=0.0                                                              1923
13990 do 13991 k=1,ni                                                      1923
      if(ixx(k).eq.0)goto 13991                                            1924
      bk=b(k,ic)                                                           1924
      gk=dot_product(r,x(:,k))                                             1925
      u=gk+xv(k,ic)*b(k,ic)                                                1925
      au=abs(u)-vp(k)*al1                                                  1926
      if(au .gt. 0.0)goto 14011                                            1926
      b(k,ic)=0.0                                                          1926
      goto 14021                                                           1927
14011 continue                                                             1928
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))   1930 
     *)
14021 continue                                                             1931
14001 continue                                                             1931
      d=b(k,ic)-bk                                                         1931
      if(abs(d).le.0.0)goto 13991                                          1932
      dlx=max(dlx,xv(k,ic)*d**2)                                           1932
      r=r-d*v*x(:,k)                                                       1933
      if(mm(k) .ne. 0)goto 14041                                           1933
      nin=nin+1                                                            1934
      if(nin .le. nx)goto 14061                                            1934
      jx=1                                                                 1934
      goto 13992                                                           1934
14061 continue                                                             1935
      mm(k)=nin                                                            1935
      m(nin)=k                                                             1936
14041 continue                                                             1937
13991 continue                                                             1938
13992 continue                                                             1938
      if(jx.gt.0)goto 13982                                                1939
      d=0.0                                                                1939
      if(intr.ne.0) d=sum(r)/xmz                                           1940
      if(d .eq. 0.0)goto 14081                                             1940
      b(0,ic)=b(0,ic)+d                                                    1940
      dlx=max(dlx,xmz*d**2)                                                1940
      r=r-d*v                                                              1940
14081 continue                                                             1941
      if(dlx.lt.shr)goto 13982                                             1942
      if(nlp .le. maxit)goto 14101                                         1942
      jerr=-ilm                                                            1942
      return                                                               1942
14101 continue                                                             1943
14110 continue                                                             1943
14111 continue                                                             1943
      nlp=nlp+1                                                            1943
      dlx=0.0                                                              1944
14120 do 14121 l=1,nin                                                     1944
      k=m(l)                                                               1944
      bk=b(k,ic)                                                           1945
      gk=dot_product(r,x(:,k))                                             1946
      u=gk+xv(k,ic)*b(k,ic)                                                1946
      au=abs(u)-vp(k)*al1                                                  1947
      if(au .gt. 0.0)goto 14141                                            1947
      b(k,ic)=0.0                                                          1947
      goto 14151                                                           1948
14141 continue                                                             1949
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))   1951 
     *)
14151 continue                                                             1952
14131 continue                                                             1952
      d=b(k,ic)-bk                                                         1952
      if(abs(d).le.0.0)goto 14121                                          1953
      dlx=max(dlx,xv(k,ic)*d**2)                                           1953
      r=r-d*v*x(:,k)                                                       1954
14121 continue                                                             1955
14122 continue                                                             1955
      d=0.0                                                                1955
      if(intr.ne.0) d=sum(r)/xmz                                           1956
      if(d .eq. 0.0)goto 14171                                             1956
      b(0,ic)=b(0,ic)+d                                                    1957
      dlx=max(dlx,xmz*d**2)                                                1957
      r=r-d*v                                                              1958
14171 continue                                                             1959
      if(dlx.lt.shr)goto 14112                                             1959
      if(nlp .le. maxit)goto 14191                                         1959
      jerr=-ilm                                                            1959
      return                                                               1959
14191 continue                                                             1960
      goto 14111                                                           1961
14112 continue                                                             1961
      goto 13981                                                           1962
13982 continue                                                             1962
      if(jx.gt.0)goto 13892                                                1963
      if(xmz*(b(0,ic)-bs(0,ic))**2.gt.shr) ix=1                            1964
      if(ix .ne. 0)goto 14211                                              1965
14220 do 14221 j=1,nin                                                     1965
      k=m(j)                                                               1966
      if(xv(k,ic)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 14241                1966
      ix=1                                                                 1966
      goto 14222                                                           1966
14241 continue                                                             1967
14221 continue                                                             1968
14222 continue                                                             1968
14211 continue                                                             1969
14250 do 14251 i=1,no                                                      1969
      fi=b(0,ic)+g(i,ic)                                                   1971
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)))         1972
      fi=min(max(exmn,fi),exmx)                                            1972
      sxp(i)=sxp(i)-q(i,ic)                                                1973
      q(i,ic)=min(max(emin*sxp(i),exp(fi)),emax*sxp(i))                    1974
      sxp(i)=sxp(i)+q(i,ic)                                                1975
14251 continue                                                             1976
14252 continue                                                             1976
13891 continue                                                             1977
13892 continue                                                             1977
      s=-sum(b(0,:))/nc                                                    1977
      b(0,:)=b(0,:)+s                                                      1977
      di=s                                                                 1978
14260 do 14261 j=1,nin                                                     1978
      l=m(j)                                                               1979
      if(vp(l) .gt. 0.0)goto 14281                                         1979
      s=sum(b(l,:))/nc                                                     1979
      goto 14291                                                           1980
14281 continue                                                             1980
      s=elc(parm,nc,cl(:,l),b(l,:),is)                                     1980
14291 continue                                                             1981
14271 continue                                                             1981
      b(l,:)=b(l,:)-s                                                      1981
      di=di-s*x(:,l)                                                       1982
14261 continue                                                             1983
14262 continue                                                             1983
      di=exp(di)                                                           1983
      sxp=sxp*di                                                           1983
14300 do 14301 ic=1,nc                                                     1983
      q(:,ic)=q(:,ic)*di                                                   1983
14301 continue                                                             1984
14302 continue                                                             1984
      if(jx.gt.0)goto 13882                                                1984
      if(ig.eq.0)goto 13882                                                1985
      if(ix .ne. 0)goto 14321                                              1986
14330 do 14331 k=1,ni                                                      1986
      if(ixx(k).eq.1)goto 14331                                            1986
      if(ju(k).eq.0)goto 14331                                             1986
      ga(k)=0.0                                                            1986
14331 continue                                                             1987
14332 continue                                                             1987
14340 do 14341 ic=1,nc                                                     1987
      r=w*(y(:,ic)-q(:,ic)/sxp)                                            1988
14350 do 14351 k=1,ni                                                      1988
      if(ixx(k).eq.1)goto 14351                                            1988
      if(ju(k).eq.0)goto 14351                                             1989
      ga(k)=max(ga(k),abs(dot_product(r,x(:,k))))                          1990
14351 continue                                                             1991
14352 continue                                                             1991
14341 continue                                                             1992
14342 continue                                                             1992
14360 do 14361 k=1,ni                                                      1992
      if(ixx(k).eq.1)goto 14361                                            1992
      if(ju(k).eq.0)goto 14361                                             1993
      if(ga(k) .le. al1*vp(k))goto 14381                                   1993
      ixx(k)=1                                                             1993
      ix=1                                                                 1993
14381 continue                                                             1994
14361 continue                                                             1995
14362 continue                                                             1995
      if(ix.eq.1) go to 11020                                              1996
      goto 13882                                                           1997
14321 continue                                                             1998
      goto 13881                                                           1999
13882 continue                                                             1999
      if(jx .le. 0)goto 14401                                              1999
      jerr=-10000-ilm                                                      1999
      goto 13802                                                           1999
14401 continue                                                             1999
      devi=0.0                                                             2000
14410 do 14411 ic=1,nc                                                     2001
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          2001
      a0(ic,ilm)=b(0,ic)                                                   2002
14420 do 14421 i=1,no                                                      2002
      if(y(i,ic).le.0.0)goto 14421                                         2003
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           2004
14421 continue                                                             2005
14422 continue                                                             2005
14411 continue                                                             2006
14412 continue                                                             2006
      kin(ilm)=nin                                                         2006
      alm(ilm)=al                                                          2006
      lmu=ilm                                                              2007
      dev(ilm)=(dev1-devi)/dev0                                            2007
      if(ig.eq.0)goto 13802                                                2008
      if(ilm.lt.mnl)goto 13801                                             2008
      if(flmin.ge.1.0)goto 13801                                           2009
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 13802             2010
      if(dev(ilm).gt.devmax)goto 13802                                     2010
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 13802                             2011
13801 continue                                                             2012
13802 continue                                                             2012
      g=log(q)                                                             2012
14430 do 14431 i=1,no                                                      2012
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         2012
14431 continue                                                             2013
14432 continue                                                             2013
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,ga,ixx)                           2014
      return                                                               2015
      end                                                                  2016
      subroutine kazero(kk,n,y,g,q,az,jerr)                                2017
      implicit double precision(a-h,o-z)                                   2018
      parameter(eps=1.0d-7)                                                2019
      double precision y(n,kk),g(n,kk),q(n),az(kk)                         2020
      double precision, dimension (:), allocatable :: s                         
      double precision, dimension (:,:), allocatable :: e                       
      allocate(e(1:n,1:kk),stat=jerr)                                           
      if(jerr.ne.0) return                                                      
      allocate(s(1:n),stat=jerr)                                           2027
      if(jerr.ne.0) return                                                 2028
      az=0.0                                                               2028
      e=exp(g)                                                             2028
14440 do 14441 i=1,n                                                       2028
      s(i)=sum(e(i,:))                                                     2028
14441 continue                                                             2029
14442 continue                                                             2029
14450 continue                                                             2029
14451 continue                                                             2029
      dm=0.0                                                               2030
14460 do 14461 k=1,kk                                                      2030
      t=0.0                                                                2030
      u=t                                                                  2031
14470 do 14471 i=1,n                                                       2031
      pik=e(i,k)/s(i)                                                      2032
      t=t+q(i)*(y(i,k)-pik)                                                2032
      u=u+q(i)*pik*(1.0-pik)                                               2033
14471 continue                                                             2034
14472 continue                                                             2034
      d=t/u                                                                2034
      az(k)=az(k)+d                                                        2034
      ed=exp(d)                                                            2034
      dm=max(dm,abs(d))                                                    2035
14480 do 14481 i=1,n                                                       2035
      z=e(i,k)                                                             2035
      e(i,k)=z*ed                                                          2035
      s(i)=s(i)-z+e(i,k)                                                   2035
14481 continue                                                             2036
14482 continue                                                             2036
14461 continue                                                             2037
14462 continue                                                             2037
      if(dm.lt.eps)goto 14452                                              2037
      goto 14451                                                           2038
14452 continue                                                             2038
      az=az-sum(az)/kk                                                     2039
      deallocate(e,s)                                                      2040
      return                                                               2041
      end                                                                  2042
      function elc(parm,n,cl,a,m)                                          2043
      implicit double precision(a-h,o-z)                                   2044
      double precision a(n),cl(2)                                          2044
      integer m(n)                                                         2045
      fn=n                                                                 2045
      am=sum(a)/fn                                                         2046
      if((parm .ne. 0.0) .and. (n .ne. 2))goto 14501                       2046
      elc=am                                                               2046
      go to 14510                                                          2046
14501 continue                                                             2047
14520 do 14521 i=1,n                                                       2047
      m(i)=i                                                               2047
14521 continue                                                             2047
14522 continue                                                             2047
      call psort7(a,m,1,n)                                                 2048
      if(a(m(1)) .ne. a(m(n)))goto 14541                                   2048
      elc=a(1)                                                             2048
      go to 14510                                                          2048
14541 continue                                                             2049
      if(mod(n,2) .ne. 1)goto 14561                                        2049
      ad=a(m(n/2+1))                                                       2049
      goto 14571                                                           2050
14561 continue                                                             2050
      ad=0.5*(a(m(n/2+1))+a(m(n/2)))                                       2050
14571 continue                                                             2051
14551 continue                                                             2051
      if(parm .ne. 1.0)goto 14591                                          2051
      elc=ad                                                               2051
      go to 14510                                                          2051
14591 continue                                                             2052
      b1=min(am,ad)                                                        2052
      b2=max(am,ad)                                                        2052
      k2=1                                                                 2053
14600 continue                                                             2053
14601 if(a(m(k2)).gt.b1)goto 14602                                         2053
      k2=k2+1                                                              2053
      goto 14601                                                           2053
14602 continue                                                             2053
      k1=k2-1                                                              2054
14610 continue                                                             2054
14611 if(a(m(k2)).ge.b2)goto 14612                                         2054
      k2=k2+1                                                              2054
      goto 14611                                                           2055
14612 continue                                                             2055
      r=parm/((1.0-parm)*fn)                                               2055
      is=0                                                                 2055
      sm=n-2*(k1-1)                                                        2056
14620 do 14621 k=k1,k2-1                                                   2056
      sm=sm-2.0                                                            2056
      s=r*sm+am                                                            2057
      if(s .le. a(m(k)) .or. s .gt. a(m(k+1)))goto 14641                   2057
      is=k                                                                 2057
      goto 14622                                                           2057
14641 continue                                                             2058
14621 continue                                                             2059
14622 continue                                                             2059
      if(is .eq. 0)goto 14661                                              2059
      elc=s                                                                2059
      go to 14510                                                          2059
14661 continue                                                             2059
      r2=2.0*r                                                             2059
      s1=a(m(k1))                                                          2059
      am2=2.0*am                                                           2060
      cri=r2*sum(abs(a-s1))+s1*(s1-am2)                                    2060
      elc=s1                                                               2061
14670 do 14671 k=k1+1,k2                                                   2061
      s=a(m(k))                                                            2061
      if(s.eq.s1)goto 14671                                                2062
      c=r2*sum(abs(a-s))+s*(s-am2)                                         2063
      if(c .ge. cri)goto 14691                                             2063
      cri=c                                                                2063
      elc=s                                                                2063
14691 continue                                                             2063
      s1=s                                                                 2064
14671 continue                                                             2065
14672 continue                                                             2065
14510 continue                                                             2065
      elc=max(maxval(a-cl(2)),min(minval(a-cl(1)),elc))                    2066
      return                                                               2067
      end                                                                  2068
      function nintot(ni,nx,nc,a,m,nin,is)                                 2069
      implicit double precision(a-h,o-z)                                   2070
      double precision a(nx,nc)                                            2070
      integer m(nx),is(ni)                                                 2071
      is=0                                                                 2071
      nintot=0                                                             2072
14700 do 14701 ic=1,nc                                                     2072
14710 do 14711 j=1,nin                                                     2072
      k=m(j)                                                               2072
      if(is(k).ne.0)goto 14711                                             2073
      if(a(j,ic).eq.0.0)goto 14711                                         2073
      is(k)=k                                                              2073
      nintot=nintot+1                                                      2074
14711 continue                                                             2074
14712 continue                                                             2074
14701 continue                                                             2075
14702 continue                                                             2075
      return                                                               2076
      end                                                                  2077
      subroutine luncomp(ni,nx,nc,ca,ia,nin,a)                             2078
      implicit double precision(a-h,o-z)                                   2079
      double precision ca(nx,nc),a(ni,nc)                                  2079
      integer ia(nx)                                                       2080
      a=0.0                                                                2081
14720 do 14721 ic=1,nc                                                     2081
      if(nin.gt.0) a(ia(1:nin),ic)=ca(1:nin,ic)                            2081
14721 continue                                                             2082
14722 continue                                                             2082
      return                                                               2083
      end                                                                  2084
      subroutine lmodval(nt,x,nc,nx,a0,ca,ia,nin,ans)                      2085
      implicit double precision(a-h,o-z)                                   2086
      double precision a0(nc),ca(nx,nc),x(nt,*),ans(nc,nt)                 2086
      integer ia(nx)                                                       2087
14730 do 14731 i=1,nt                                                      2087
14740 do 14741 ic=1,nc                                                     2087
      ans(ic,i)=a0(ic)                                                     2089
      if(nin.gt.0) ans(ic,i)=ans(ic,i)+dot_product(ca(1:nin,ic),x(i,ia(1   2090 
     *:nin)))
14741 continue                                                             2090
14742 continue                                                             2090
14731 continue                                                             2091
14732 continue                                                             2091
      return                                                               2092
      end                                                                  2093
      subroutine splognet (parm,no,ni,nc,x,ix,jx,y,g,jd,vp,cl,ne,nx,nlam   2095 
     *,flmin,  ulam,thr,isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,al
     *m,nlp,jerr)
      implicit double precision(a-h,o-z)                                   2096
      double precision x(*),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nlam)     2097
      double precision ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl   2098 
     *(2,ni)
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           2099
      double precision, dimension (:), allocatable :: xm,xs,ww,vq,xv            
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 14761                                    2103
      jerr=10000                                                           2103
      return                                                               2103
14761 continue                                                             2104
      allocate(ww(1:no),stat=jerr)                                         2105
      if(jerr.ne.0) return                                                 2106
      allocate(ju(1:ni),stat=jerr)                                         2107
      if(jerr.ne.0) return                                                 2108
      allocate(vq(1:ni),stat=jerr)                                         2109
      if(jerr.ne.0) return                                                 2110
      allocate(xm(1:ni),stat=jerr)                                         2111
      if(jerr.ne.0) return                                                 2112
      allocate(xs(1:ni),stat=jerr)                                         2113
      if(jerr.ne.0) return                                                 2114
      if(kopt .ne. 2)goto 14781                                            2114
      allocate(xv(1:ni),stat=jerr)                                         2114
      if(jerr.ne.0) return                                                 2114
14781 continue                                                             2116
      call spchkvars(no,ni,x,ix,ju)                                        2117
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2118
      if(maxval(ju) .gt. 0)goto 14801                                      2118
      jerr=7777                                                            2118
      return                                                               2118
14801 continue                                                             2119
      vq=max(0d0,vp)                                                       2119
      vq=vq*ni/sum(vq)                                                     2120
14810 do 14811 i=1,no                                                      2120
      ww(i)=sum(y(i,:))                                                    2120
      if(ww(i).gt.0.0) y(i,:)=y(i,:)/ww(i)                                 2120
14811 continue                                                             2121
14812 continue                                                             2121
      sw=sum(ww)                                                           2121
      ww=ww/sw                                                             2122
      if(nc .ne. 1)goto 14831                                              2122
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs)                2123
      if(isd .le. 0)goto 14851                                             2123
14860 do 14861 j=1,ni                                                      2123
      cl(:,j)=cl(:,j)*xs(j)                                                2123
14861 continue                                                             2123
14862 continue                                                             2123
14851 continue                                                             2124
      call sprlognet2n(parm,no,ni,x,ix,jx,y(:,1),g(:,1),ww,ju,vq,cl,ne,n   2127 
     *x,nlam,  flmin,ulam,thr,isd,intr,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin
     *,dev0,dev,  alm,nlp,jerr)
      goto 14821                                                           2128
14831 if(kopt .ne. 2)goto 14871                                            2129
      call multsplstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs,xv)         2130
      if(isd .le. 0)goto 14891                                             2130
14900 do 14901 j=1,ni                                                      2130
      cl(:,j)=cl(:,j)*xs(j)                                                2130
14901 continue                                                             2130
14902 continue                                                             2130
14891 continue                                                             2131
      call multsprlognetn(parm,no,ni,nc,x,ix,jx,y,g,ww,ju,vq,cl,ne,nx,nl   2133 
     *am,flmin,  ulam,thr,intr,maxit,xv,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,
     *alm,nlp,jerr)
      goto 14911                                                           2134
14871 continue                                                             2134
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs)                2135
      if(isd .le. 0)goto 14931                                             2135
14940 do 14941 j=1,ni                                                      2135
      cl(:,j)=cl(:,j)*xs(j)                                                2135
14941 continue                                                             2135
14942 continue                                                             2135
14931 continue                                                             2136
      call sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,ww,ju,vq,cl,ne,nx,nlam,f   2139 
     *lmin,  ulam,thr,isd,intr,maxit,kopt,xm,xs,lmu,a0,ca,  ia,nin,dev0,
     *dev,alm,nlp,jerr)
14911 continue                                                             2140
14821 continue                                                             2140
      if(jerr.gt.0) return                                                 2140
      dev0=2.0*sw*dev0                                                     2141
14950 do 14951 k=1,lmu                                                     2141
      nk=nin(k)                                                            2142
14960 do 14961 ic=1,nc                                                     2142
      if(isd .le. 0)goto 14981                                             2142
14990 do 14991 l=1,nk                                                      2142
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                      2142
14991 continue                                                             2142
14992 continue                                                             2142
14981 continue                                                             2143
      if(intr .ne. 0)goto 15011                                            2143
      a0(ic,k)=0.0                                                         2143
      goto 15021                                                           2144
15011 continue                                                             2144
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))            2144
15021 continue                                                             2145
15001 continue                                                             2145
14961 continue                                                             2146
14962 continue                                                             2146
14951 continue                                                             2147
14952 continue                                                             2147
      deallocate(ww,ju,vq,xm,xs)                                           2147
      if(kopt.eq.2) deallocate(xv)                                         2148
      return                                                               2149
      end                                                                  2150
      subroutine multsplstandard2(no,ni,x,ix,jx,w,ju,isd,intr,xm,xs,xv)    2151
      implicit double precision(a-h,o-z)                                   2152
      double precision x(*),w(no),xm(ni),xs(ni),xv(ni)                     2152
      integer ix(*),jx(*),ju(ni)                                           2153
      if(intr .ne. 0)goto 15041                                            2154
15050 do 15051 j=1,ni                                                      2154
      if(ju(j).eq.0)goto 15051                                             2154
      xm(j)=0.0                                                            2154
      jb=ix(j)                                                             2154
      je=ix(j+1)-1                                                         2155
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)                          2156
      if(isd .eq. 0)goto 15071                                             2156
      xbq=dot_product(w(jx(jb:je)),x(jb:je))**2                            2156
      vc=xv(j)-xbq                                                         2157
      xs(j)=sqrt(vc)                                                       2157
      xv(j)=1.0+xbq/vc                                                     2158
      goto 15081                                                           2159
15071 continue                                                             2159
      xs(j)=1.0                                                            2159
15081 continue                                                             2160
15061 continue                                                             2160
15051 continue                                                             2161
15052 continue                                                             2161
      return                                                               2162
15041 continue                                                             2163
15090 do 15091 j=1,ni                                                      2163
      if(ju(j).eq.0)goto 15091                                             2163
      jb=ix(j)                                                             2163
      je=ix(j+1)-1                                                         2164
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2165
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 2166
      if(isd .le. 0)goto 15111                                             2166
      xs(j)=sqrt(xv(j))                                                    2166
      xv(j)=1.0                                                            2166
15111 continue                                                             2167
15091 continue                                                             2168
15092 continue                                                             2168
      if(isd.eq.0) xs=1.0                                                  2169
      return                                                               2170
      end                                                                  2171
      subroutine splstandard2(no,ni,x,ix,jx,w,ju,isd,intr,xm,xs)           2172
      implicit double precision(a-h,o-z)                                   2173
      double precision x(*),w(no),xm(ni),xs(ni)                            2173
      integer ix(*),jx(*),ju(ni)                                           2174
      if(intr .ne. 0)goto 15131                                            2175
15140 do 15141 j=1,ni                                                      2175
      if(ju(j).eq.0)goto 15141                                             2175
      xm(j)=0.0                                                            2175
      jb=ix(j)                                                             2175
      je=ix(j+1)-1                                                         2176
      if(isd .eq. 0)goto 15161                                             2177
      vc=dot_product(w(jx(jb:je)),x(jb:je)**2)  -dot_product(w(jx(jb:je)   2179 
     *),x(jb:je))**2
      xs(j)=sqrt(vc)                                                       2180
      goto 15171                                                           2181
15161 continue                                                             2181
      xs(j)=1.0                                                            2181
15171 continue                                                             2182
15151 continue                                                             2182
15141 continue                                                             2183
15142 continue                                                             2183
      return                                                               2184
15131 continue                                                             2185
15180 do 15181 j=1,ni                                                      2185
      if(ju(j).eq.0)goto 15181                                             2185
      jb=ix(j)                                                             2185
      je=ix(j+1)-1                                                         2186
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2187
      if(isd.ne.0) xs(j)=sqrt(dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j   2188 
     *)**2)
15181 continue                                                             2189
15182 continue                                                             2189
      if(isd.eq.0) xs=1.0                                                  2190
      return                                                               2191
      end                                                                  2192
      subroutine sprlognet2n (parm,no,ni,x,ix,jx,y,g,w,ju,vp,cl,ne,nx,nl   2195 
     *am,  flmin,ulam,shri,isd,intr,maxit,kopt,xb,xs,  lmu,a0,a,m,kin,de
     *v0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   2196
      double precision x(*),y(no),g(no),w(no),vp(ni),ulam(nlam),cl(2,ni)   2197
      double precision a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)             2198
      double precision xb(ni),xs(ni)                                       2198
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           2199
      double precision, dimension (:), allocatable :: xm,b,bs,v,r               
      double precision, dimension (:), allocatable :: sc,xv,q,ga                
      integer, dimension (:), allocatable :: mm,ixx                             
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               2205
      allocate(b(0:ni),stat=jerr)                                          2206
      if(jerr.ne.0) return                                                 2207
      allocate(xm(0:ni),stat=jerr)                                         2208
      if(jerr.ne.0) return                                                 2209
      allocate(xv(1:ni),stat=jerr)                                         2210
      if(jerr.ne.0) return                                                 2211
      allocate(bs(0:ni),stat=jerr)                                         2212
      if(jerr.ne.0) return                                                 2213
      allocate(ga(1:ni),stat=jerr)                                         2214
      if(jerr.ne.0) return                                                 2215
      allocate(mm(1:ni),stat=jerr)                                         2216
      if(jerr.ne.0) return                                                 2217
      allocate(ixx(1:ni),stat=jerr)                                        2218
      if(jerr.ne.0) return                                                 2219
      allocate(q(1:no),stat=jerr)                                          2220
      if(jerr.ne.0) return                                                 2221
      allocate(r(1:no),stat=jerr)                                          2222
      if(jerr.ne.0) return                                                 2223
      allocate(v(1:no),stat=jerr)                                          2224
      if(jerr.ne.0) return                                                 2225
      allocate(sc(1:no),stat=jerr)                                         2226
      if(jerr.ne.0) return                                                 2227
      fmax=log(1.0/pmin-1.0)                                               2227
      fmin=-fmax                                                           2227
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                      2228
      bta=parm                                                             2228
      omb=1.0-bta                                                          2229
      q0=dot_product(w,y)                                                  2229
      if(q0 .gt. pmin)goto 15201                                           2229
      jerr=8001                                                            2229
      return                                                               2229
15201 continue                                                             2230
      if(q0 .lt. 1.0-pmin)goto 15221                                       2230
      jerr=9001                                                            2230
      return                                                               2230
15221 continue                                                             2231
      if(intr.eq.0) q0=0.5                                                 2231
      bz=0.0                                                               2231
      if(intr.ne.0) bz=log(q0/(1.0-q0))                                    2232
      if(nonzero(no,g) .ne. 0)goto 15241                                   2232
      vi=q0*(1.0-q0)                                                       2232
      b(0)=bz                                                              2232
      v=vi*w                                                               2233
      r=w*(y-q0)                                                           2233
      q=q0                                                                 2233
      xm(0)=vi                                                             2233
      dev1=-(bz*q0+log(1.0-q0))                                            2234
      goto 15251                                                           2235
15241 continue                                                             2235
      b(0)=0.0                                                             2236
      if(intr .eq. 0)goto 15271                                            2236
      b(0)=azero(no,y,g,w,jerr)                                            2236
      if(jerr.ne.0) return                                                 2236
15271 continue                                                             2237
      q=1.0/(1.0+exp(-b(0)-g))                                             2237
      v=w*q*(1.0-q)                                                        2237
      r=w*(y-q)                                                            2237
      xm(0)=sum(v)                                                         2238
      dev1=-(b(0)*q0+dot_product(w,y*g+log(1.0-q)))                        2239
15251 continue                                                             2240
15231 continue                                                             2240
      if(kopt .le. 0)goto 15291                                            2241
      if(isd .le. 0 .or. intr .eq. 0)goto 15311                            2241
      xv=0.25                                                              2241
      goto 15321                                                           2242
15311 continue                                                             2243
15330 do 15331 j=1,ni                                                      2243
      if(ju(j).eq.0)goto 15331                                             2243
      jb=ix(j)                                                             2243
      je=ix(j+1)-1                                                         2244
      xv(j)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)          2245
15331 continue                                                             2246
15332 continue                                                             2246
15321 continue                                                             2247
15301 continue                                                             2247
15291 continue                                                             2248
      b(1:ni)=0.0                                                          2248
      dev0=dev1                                                            2249
15340 do 15341 i=1,no                                                      2249
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                        2250
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))              2251
15341 continue                                                             2253
15342 continue                                                             2253
      alf=1.0                                                              2255
      if(flmin .ge. 1.0)goto 15361                                         2255
      eqs=max(eps,flmin)                                                   2255
      alf=eqs**(1.0/(nlam-1))                                              2255
15361 continue                                                             2256
      m=0                                                                  2256
      mm=0                                                                 2256
      nin=0                                                                2256
      o=0.0                                                                2256
      svr=o                                                                2256
      mnl=min(mnlam,nlam)                                                  2256
      bs=0.0                                                               2256
      nlp=0                                                                2256
      nin=nlp                                                              2257
      shr=shri*dev0                                                        2257
      al=0.0                                                               2257
      ixx=0                                                                2258
15370 do 15371 j=1,ni                                                      2258
      if(ju(j).eq.0)goto 15371                                             2259
      jb=ix(j)                                                             2259
      je=ix(j+1)-1                                                         2259
      jn=ix(j+1)-ix(j)                                                     2260
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 2261
      gj=dot_product(sc(1:jn),x(jb:je))                                    2262
      ga(j)=abs((gj-svr*xb(j))/xs(j))                                      2263
15371 continue                                                             2264
15372 continue                                                             2264
15380 do 15381 ilm=1,nlam                                                  2264
      al0=al                                                               2265
      if(flmin .lt. 1.0)goto 15401                                         2265
      al=ulam(ilm)                                                         2265
      goto 15391                                                           2266
15401 if(ilm .le. 2)goto 15411                                             2266
      al=al*alf                                                            2266
      goto 15391                                                           2267
15411 if(ilm .ne. 1)goto 15421                                             2267
      al=big                                                               2267
      goto 15431                                                           2268
15421 continue                                                             2268
      al0=0.0                                                              2269
15440 do 15441 j=1,ni                                                      2269
      if(ju(j).eq.0)goto 15441                                             2269
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2269
15441 continue                                                             2270
15442 continue                                                             2270
      al0=al0/max(bta,1.0d-3)                                              2270
      al=alf*al0                                                           2271
15431 continue                                                             2272
15391 continue                                                             2272
      al2=al*omb                                                           2272
      al1=al*bta                                                           2272
      tlam=bta*(2.0*al-al0)                                                2273
15450 do 15451 k=1,ni                                                      2273
      if(ixx(k).eq.1)goto 15451                                            2273
      if(ju(k).eq.0)goto 15451                                             2274
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     2275
15451 continue                                                             2276
15452 continue                                                             2276
11020 continue                                                             2277
15460 continue                                                             2277
15461 continue                                                             2277
      bs(0)=b(0)                                                           2277
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                2278
15470 do 15471 j=1,ni                                                      2278
      if(ixx(j).eq.0)goto 15471                                            2279
      jb=ix(j)                                                             2279
      je=ix(j+1)-1                                                         2279
      jn=ix(j+1)-ix(j)                                                     2280
      sc(1:jn)=v(jx(jb:je))                                                2281
      xm(j)=dot_product(sc(1:jn),x(jb:je))                                 2282
      if(kopt .ne. 0)goto 15491                                            2283
      xv(j)=dot_product(sc(1:jn),x(jb:je)**2)                              2284
      xv(j)=(xv(j)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2                2285
15491 continue                                                             2286
15471 continue                                                             2287
15472 continue                                                             2287
15500 continue                                                             2287
15501 continue                                                             2287
      nlp=nlp+1                                                            2287
      dlx=0.0                                                              2288
15510 do 15511 k=1,ni                                                      2288
      if(ixx(k).eq.0)goto 15511                                            2289
      jb=ix(k)                                                             2289
      je=ix(k+1)-1                                                         2289
      jn=ix(k+1)-ix(k)                                                     2289
      bk=b(k)                                                              2290
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 2291
      gk=dot_product(sc(1:jn),x(jb:je))                                    2292
      gk=(gk-svr*xb(k))/xs(k)                                              2293
      u=gk+xv(k)*b(k)                                                      2293
      au=abs(u)-vp(k)*al1                                                  2294
      if(au .gt. 0.0)goto 15531                                            2294
      b(k)=0.0                                                             2294
      goto 15541                                                           2295
15531 continue                                                             2296
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))          2297
15541 continue                                                             2298
15521 continue                                                             2298
      d=b(k)-bk                                                            2298
      if(abs(d).le.0.0)goto 15511                                          2298
      dlx=max(dlx,xv(k)*d**2)                                              2299
      if(mm(k) .ne. 0)goto 15561                                           2299
      nin=nin+1                                                            2299
      if(nin.gt.nx)goto 15512                                              2300
      mm(k)=nin                                                            2300
      m(nin)=k                                                             2300
      sc(1:jn)=v(jx(jb:je))                                                2301
      xm(k)=dot_product(sc(1:jn),x(jb:je))                                 2302
15561 continue                                                             2303
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              2304
      o=o+d*(xb(k)/xs(k))                                                  2305
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  2306
15511 continue                                                             2307
15512 continue                                                             2307
      if(nin.gt.nx)goto 15502                                              2308
      d=0.0                                                                2308
      if(intr.ne.0) d=svr/xm(0)                                            2309
      if(d .eq. 0.0)goto 15581                                             2309
      b(0)=b(0)+d                                                          2309
      dlx=max(dlx,xm(0)*d**2)                                              2309
      r=r-d*v                                                              2310
      svr=svr-d*xm(0)                                                      2311
15581 continue                                                             2312
      if(dlx.lt.shr)goto 15502                                             2313
      if(nlp .le. maxit)goto 15601                                         2313
      jerr=-ilm                                                            2313
      return                                                               2313
15601 continue                                                             2314
15610 continue                                                             2314
15611 continue                                                             2314
      nlp=nlp+1                                                            2314
      dlx=0.0                                                              2315
15620 do 15621 l=1,nin                                                     2315
      k=m(l)                                                               2315
      jb=ix(k)                                                             2315
      je=ix(k+1)-1                                                         2316
      jn=ix(k+1)-ix(k)                                                     2316
      bk=b(k)                                                              2317
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 2318
      gk=dot_product(sc(1:jn),x(jb:je))                                    2319
      gk=(gk-svr*xb(k))/xs(k)                                              2320
      u=gk+xv(k)*b(k)                                                      2320
      au=abs(u)-vp(k)*al1                                                  2321
      if(au .gt. 0.0)goto 15641                                            2321
      b(k)=0.0                                                             2321
      goto 15651                                                           2322
15641 continue                                                             2323
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))          2324
15651 continue                                                             2325
15631 continue                                                             2325
      d=b(k)-bk                                                            2325
      if(abs(d).le.0.0)goto 15621                                          2325
      dlx=max(dlx,xv(k)*d**2)                                              2326
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              2327
      o=o+d*(xb(k)/xs(k))                                                  2328
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  2329
15621 continue                                                             2330
15622 continue                                                             2330
      d=0.0                                                                2330
      if(intr.ne.0) d=svr/xm(0)                                            2331
      if(d .eq. 0.0)goto 15671                                             2331
      b(0)=b(0)+d                                                          2331
      dlx=max(dlx,xm(0)*d**2)                                              2331
      r=r-d*v                                                              2332
      svr=svr-d*xm(0)                                                      2333
15671 continue                                                             2334
      if(dlx.lt.shr)goto 15612                                             2335
      if(nlp .le. maxit)goto 15691                                         2335
      jerr=-ilm                                                            2335
      return                                                               2335
15691 continue                                                             2336
      goto 15611                                                           2337
15612 continue                                                             2337
      goto 15501                                                           2338
15502 continue                                                             2338
      if(nin.gt.nx)goto 15462                                              2339
      sc=b(0)                                                              2339
      b0=0.0                                                               2340
15700 do 15701 j=1,nin                                                     2340
      l=m(j)                                                               2340
      jb=ix(l)                                                             2340
      je=ix(l+1)-1                                                         2341
      sc(jx(jb:je))=sc(jx(jb:je))+b(l)*x(jb:je)/xs(l)                      2342
      b0=b0-b(l)*xb(l)/xs(l)                                               2343
15701 continue                                                             2344
15702 continue                                                             2344
      sc=sc+b0                                                             2345
15710 do 15711 i=1,no                                                      2345
      fi=sc(i)+g(i)                                                        2346
      if(fi .ge. fmin)goto 15731                                           2346
      q(i)=0.0                                                             2346
      goto 15721                                                           2346
15731 if(fi .le. fmax)goto 15741                                           2346
      q(i)=1.0                                                             2346
      goto 15751                                                           2347
15741 continue                                                             2347
      q(i)=1.0/(1.0+exp(-fi))                                              2347
15751 continue                                                             2348
15721 continue                                                             2348
15711 continue                                                             2349
15712 continue                                                             2349
      v=w*q*(1.0-q)                                                        2349
      xm(0)=sum(v)                                                         2349
      if(xm(0).lt.vmin)goto 15462                                          2350
      r=w*(y-q)                                                            2350
      svr=sum(r)                                                           2350
      o=0.0                                                                2351
      if(xm(0)*(b(0)-bs(0))**2 .ge. shr)goto 15771                         2351
      kx=0                                                                 2352
15780 do 15781 j=1,nin                                                     2352
      k=m(j)                                                               2353
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 15781                           2353
      kx=1                                                                 2353
      goto 15782                                                           2354
15781 continue                                                             2355
15782 continue                                                             2355
      if(kx .ne. 0)goto 15801                                              2356
15810 do 15811 j=1,ni                                                      2356
      if(ixx(j).eq.1)goto 15811                                            2356
      if(ju(j).eq.0)goto 15811                                             2357
      jb=ix(j)                                                             2357
      je=ix(j+1)-1                                                         2357
      jn=ix(j+1)-ix(j)                                                     2358
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 2359
      gj=dot_product(sc(1:jn),x(jb:je))                                    2360
      ga(j)=abs((gj-svr*xb(j))/xs(j))                                      2361
      if(ga(j) .le. al1*vp(j))goto 15831                                   2361
      ixx(j)=1                                                             2361
      kx=1                                                                 2361
15831 continue                                                             2362
15811 continue                                                             2363
15812 continue                                                             2363
      if(kx.eq.1) go to 11020                                              2364
      goto 15462                                                           2365
15801 continue                                                             2366
15771 continue                                                             2367
      goto 15461                                                           2368
15462 continue                                                             2368
      if(nin .le. nx)goto 15851                                            2368
      jerr=-10000-ilm                                                      2368
      goto 15382                                                           2368
15851 continue                                                             2369
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                2369
      kin(ilm)=nin                                                         2370
      a0(ilm)=b(0)                                                         2370
      alm(ilm)=al                                                          2370
      lmu=ilm                                                              2371
      devi=dev2(no,w,y,q,pmin)                                             2372
      dev(ilm)=(dev1-devi)/dev0                                            2373
      if(ilm.lt.mnl)goto 15381                                             2373
      if(flmin.ge.1.0)goto 15381                                           2374
      me=0                                                                 2374
15860 do 15861 j=1,nin                                                     2374
      if(a(j,ilm).ne.0.0) me=me+1                                          2374
15861 continue                                                             2374
15862 continue                                                             2374
      if(me.gt.ne)goto 15382                                               2375
      if(dev(ilm).gt.devmax)goto 15382                                     2375
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 15382                             2376
      if(xm(0).lt.vmin)goto 15382                                          2377
15381 continue                                                             2378
15382 continue                                                             2378
      g=log(q/(1.0-q))                                                     2379
      deallocate(xm,b,bs,v,r,sc,xv,q,mm,ga,ixx)                            2380
      return                                                               2381
      end                                                                  2382
      subroutine sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,w,ju,vp,cl,ne,nx,n   2384 
     *lam,flmin,  ulam,shri,isd,intr,maxit,kopt,xb,xs,lmu,a0,a,m,kin,dev
     *0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   2385
      double precision x(*),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam),xb   2386 
     *(ni),xs(ni)
      double precision a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(   2387 
     *2,ni)
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           2388
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
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               2403
      exmn=-exmx                                                           2404
      allocate(xm(0:ni),stat=jerr)                                         2405
      if(jerr.ne.0) return                                                 2406
      allocate(r(1:no),stat=jerr)                                          2407
      if(jerr.ne.0) return                                                 2408
      allocate(v(1:no),stat=jerr)                                          2409
      if(jerr.ne.0) return                                                 2410
      allocate(mm(1:ni),stat=jerr)                                         2411
      if(jerr.ne.0) return                                                 2412
      allocate(ga(1:ni),stat=jerr)                                         2413
      if(jerr.ne.0) return                                                 2414
      allocate(iy(1:ni),stat=jerr)                                         2415
      if(jerr.ne.0) return                                                 2416
      allocate(is(1:max(nc,ni)),stat=jerr)                                 2417
      if(jerr.ne.0) return                                                 2418
      allocate(sxp(1:no),stat=jerr)                                        2419
      if(jerr.ne.0) return                                                 2420
      allocate(sxpl(1:no),stat=jerr)                                       2421
      if(jerr.ne.0) return                                                 2422
      allocate(sc(1:no),stat=jerr)                                         2423
      if(jerr.ne.0) return                                                 2424
      pmax=1.0-pmin                                                        2424
      emin=pmin/pmax                                                       2424
      emax=1.0/emin                                                        2425
      pfm=(1.0+pmin)*pmin                                                  2425
      pfx=(1.0-pmin)*pmax                                                  2425
      vmin=pfm*pmax                                                        2426
      bta=parm                                                             2426
      omb=1.0-bta                                                          2426
      dev1=0.0                                                             2426
      dev0=0.0                                                             2427
15870 do 15871 ic=1,nc                                                     2427
      q0=dot_product(w,y(:,ic))                                            2428
      if(q0 .gt. pmin)goto 15891                                           2428
      jerr =8000+ic                                                        2428
      return                                                               2428
15891 continue                                                             2429
      if(q0 .lt. 1.0-pmin)goto 15911                                       2429
      jerr =9000+ic                                                        2429
      return                                                               2429
15911 continue                                                             2430
      if(intr.eq.0) q0=1.0/nc                                              2431
      b(1:ni,ic)=0.0                                                       2431
      b(0,ic)=0.0                                                          2432
      if(intr .eq. 0)goto 15931                                            2432
      b(0,ic)=log(q0)                                                      2432
      dev1=dev1-q0*b(0,ic)                                                 2432
15931 continue                                                             2433
15871 continue                                                             2434
15872 continue                                                             2434
      if(intr.eq.0) dev1=log(float(nc))                                    2434
      iy=0                                                                 2434
      al=0.0                                                               2435
      if(nonzero(no*nc,g) .ne. 0)goto 15951                                2436
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         2436
      sxp=0.0                                                              2437
15960 do 15961 ic=1,nc                                                     2437
      q(:,ic)=exp(b(0,ic))                                                 2437
      sxp=sxp+q(:,ic)                                                      2437
15961 continue                                                             2438
15962 continue                                                             2438
      goto 15971                                                           2439
15951 continue                                                             2439
15980 do 15981 i=1,no                                                      2439
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         2439
15981 continue                                                             2439
15982 continue                                                             2439
      sxp=0.0                                                              2440
      if(intr .ne. 0)goto 16001                                            2440
      b(0,:)=0.0                                                           2440
      goto 16011                                                           2441
16001 continue                                                             2441
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 2441
      if(jerr.ne.0) return                                                 2441
16011 continue                                                             2442
15991 continue                                                             2442
      dev1=0.0                                                             2443
16020 do 16021 ic=1,nc                                                     2443
      q(:,ic)=b(0,ic)+g(:,ic)                                              2444
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             2445
      q(:,ic)=exp(q(:,ic))                                                 2445
      sxp=sxp+q(:,ic)                                                      2446
16021 continue                                                             2447
16022 continue                                                             2447
      sxpl=w*log(sxp)                                                      2447
16030 do 16031 ic=1,nc                                                     2447
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  2447
16031 continue                                                             2448
16032 continue                                                             2448
15971 continue                                                             2449
15941 continue                                                             2449
16040 do 16041 ic=1,nc                                                     2449
16050 do 16051 i=1,no                                                      2449
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               2449
16051 continue                                                             2449
16052 continue                                                             2449
16041 continue                                                             2450
16042 continue                                                             2450
      dev0=dev0+dev1                                                       2451
      if(kopt .le. 0)goto 16071                                            2452
      if(isd .le. 0 .or. intr .eq. 0)goto 16091                            2452
      xv=0.25                                                              2452
      goto 16101                                                           2453
16091 continue                                                             2454
16110 do 16111 j=1,ni                                                      2454
      if(ju(j).eq.0)goto 16111                                             2454
      jb=ix(j)                                                             2454
      je=ix(j+1)-1                                                         2455
      xv(j,:)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)        2456
16111 continue                                                             2457
16112 continue                                                             2457
16101 continue                                                             2458
16081 continue                                                             2458
16071 continue                                                             2460
      alf=1.0                                                              2462
      if(flmin .ge. 1.0)goto 16131                                         2462
      eqs=max(eps,flmin)                                                   2462
      alf=eqs**(1.0/(nlam-1))                                              2462
16131 continue                                                             2463
      m=0                                                                  2463
      mm=0                                                                 2463
      nin=0                                                                2463
      nlp=0                                                                2463
      mnl=min(mnlam,nlam)                                                  2463
      bs=0.0                                                               2463
      svr=0.0                                                              2463
      o=0.0                                                                2464
      shr=shri*dev0                                                        2464
      ga=0.0                                                               2465
16140 do 16141 ic=1,nc                                                     2465
      v=q(:,ic)/sxp                                                        2465
      r=w*(y(:,ic)-v)                                                      2465
      v=w*v*(1.0-v)                                                        2466
16150 do 16151 j=1,ni                                                      2466
      if(ju(j).eq.0)goto 16151                                             2467
      jb=ix(j)                                                             2467
      je=ix(j+1)-1                                                         2467
      jn=ix(j+1)-ix(j)                                                     2468
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2469
      gj=dot_product(sc(1:jn),x(jb:je))                                    2470
      ga(j)=max(ga(j),abs(gj-svr*xb(j))/xs(j))                             2471
16151 continue                                                             2472
16152 continue                                                             2472
16141 continue                                                             2473
16142 continue                                                             2473
16160 do 16161 ilm=1,nlam                                                  2473
      al0=al                                                               2474
      if(flmin .lt. 1.0)goto 16181                                         2474
      al=ulam(ilm)                                                         2474
      goto 16171                                                           2475
16181 if(ilm .le. 2)goto 16191                                             2475
      al=al*alf                                                            2475
      goto 16171                                                           2476
16191 if(ilm .ne. 1)goto 16201                                             2476
      al=big                                                               2476
      goto 16211                                                           2477
16201 continue                                                             2477
      al0=0.0                                                              2478
16220 do 16221 j=1,ni                                                      2478
      if(ju(j).eq.0)goto 16221                                             2478
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2478
16221 continue                                                             2479
16222 continue                                                             2479
      al0=al0/max(bta,1.0d-3)                                              2479
      al=alf*al0                                                           2480
16211 continue                                                             2481
16171 continue                                                             2481
      al2=al*omb                                                           2481
      al1=al*bta                                                           2481
      tlam=bta*(2.0*al-al0)                                                2482
16230 do 16231 k=1,ni                                                      2482
      if(iy(k).eq.1)goto 16231                                             2482
      if(ju(k).eq.0)goto 16231                                             2483
      if(ga(k).gt.tlam*vp(k)) iy(k)=1                                      2484
16231 continue                                                             2485
16232 continue                                                             2485
11020 continue                                                             2486
16240 continue                                                             2486
16241 continue                                                             2486
      ixx=0                                                                2486
      jxx=ixx                                                              2486
      ig=0                                                                 2487
16250 do 16251 ic=1,nc                                                     2487
      bs(0,ic)=b(0,ic)                                                     2488
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          2489
      xm(0)=0.0                                                            2489
      svr=0.0                                                              2489
      o=0.0                                                                2490
16260 do 16261 i=1,no                                                      2490
      pic=q(i,ic)/sxp(i)                                                   2491
      if(pic .ge. pfm)goto 16281                                           2491
      pic=0.0                                                              2491
      v(i)=0.0                                                             2491
      goto 16271                                                           2492
16281 if(pic .le. pfx)goto 16291                                           2492
      pic=1.0                                                              2492
      v(i)=0.0                                                             2492
      goto 16301                                                           2493
16291 continue                                                             2493
      v(i)=w(i)*pic*(1.0-pic)                                              2493
      xm(0)=xm(0)+v(i)                                                     2493
16301 continue                                                             2494
16271 continue                                                             2494
      r(i)=w(i)*(y(i,ic)-pic)                                              2494
      svr=svr+r(i)                                                         2495
16261 continue                                                             2496
16262 continue                                                             2496
      if(xm(0).le.vmin)goto 16251                                          2496
      ig=1                                                                 2497
16310 do 16311 j=1,ni                                                      2497
      if(iy(j).eq.0)goto 16311                                             2498
      jb=ix(j)                                                             2498
      je=ix(j+1)-1                                                         2499
      xm(j)=dot_product(v(jx(jb:je)),x(jb:je))                             2500
      if(kopt .ne. 0)goto 16331                                            2501
      xv(j,ic)=dot_product(v(jx(jb:je)),x(jb:je)**2)                       2502
      xv(j,ic)=(xv(j,ic)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2          2503
16331 continue                                                             2504
16311 continue                                                             2505
16312 continue                                                             2505
16340 continue                                                             2505
16341 continue                                                             2505
      nlp=nlp+1                                                            2505
      dlx=0.0                                                              2506
16350 do 16351 k=1,ni                                                      2506
      if(iy(k).eq.0)goto 16351                                             2507
      jb=ix(k)                                                             2507
      je=ix(k+1)-1                                                         2507
      jn=ix(k+1)-ix(k)                                                     2507
      bk=b(k,ic)                                                           2508
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2509
      gk=dot_product(sc(1:jn),x(jb:je))                                    2510
      gk=(gk-svr*xb(k))/xs(k)                                              2511
      u=gk+xv(k,ic)*b(k,ic)                                                2511
      au=abs(u)-vp(k)*al1                                                  2512
      if(au .gt. 0.0)goto 16371                                            2512
      b(k,ic)=0.0                                                          2512
      goto 16381                                                           2513
16371 continue                                                             2514
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))   2516 
     *)
16381 continue                                                             2517
16361 continue                                                             2517
      d=b(k,ic)-bk                                                         2517
      if(abs(d).le.0.0)goto 16351                                          2518
      dlx=max(dlx,xv(k,ic)*d**2)                                           2519
      if(mm(k) .ne. 0)goto 16401                                           2519
      nin=nin+1                                                            2520
      if(nin .le. nx)goto 16421                                            2520
      jxx=1                                                                2520
      goto 16352                                                           2520
16421 continue                                                             2521
      mm(k)=nin                                                            2521
      m(nin)=k                                                             2522
      xm(k)=dot_product(v(jx(jb:je)),x(jb:je))                             2523
16401 continue                                                             2524
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              2525
      o=o+d*(xb(k)/xs(k))                                                  2526
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  2527
16351 continue                                                             2528
16352 continue                                                             2528
      if(jxx.gt.0)goto 16342                                               2529
      d=0.0                                                                2529
      if(intr.ne.0) d=svr/xm(0)                                            2530
      if(d .eq. 0.0)goto 16441                                             2530
      b(0,ic)=b(0,ic)+d                                                    2530
      dlx=max(dlx,xm(0)*d**2)                                              2531
      r=r-d*v                                                              2531
      svr=svr-d*xm(0)                                                      2532
16441 continue                                                             2533
      if(dlx.lt.shr)goto 16342                                             2533
      if(nlp .le. maxit)goto 16461                                         2533
      jerr=-ilm                                                            2533
      return                                                               2533
16461 continue                                                             2534
16470 continue                                                             2534
16471 continue                                                             2534
      nlp=nlp+1                                                            2534
      dlx=0.0                                                              2535
16480 do 16481 l=1,nin                                                     2535
      k=m(l)                                                               2535
      jb=ix(k)                                                             2535
      je=ix(k+1)-1                                                         2536
      jn=ix(k+1)-ix(k)                                                     2536
      bk=b(k,ic)                                                           2537
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2538
      gk=dot_product(sc(1:jn),x(jb:je))                                    2539
      gk=(gk-svr*xb(k))/xs(k)                                              2540
      u=gk+xv(k,ic)*b(k,ic)                                                2540
      au=abs(u)-vp(k)*al1                                                  2541
      if(au .gt. 0.0)goto 16501                                            2541
      b(k,ic)=0.0                                                          2541
      goto 16511                                                           2542
16501 continue                                                             2543
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))   2545 
     *)
16511 continue                                                             2546
16491 continue                                                             2546
      d=b(k,ic)-bk                                                         2546
      if(abs(d).le.0.0)goto 16481                                          2547
      dlx=max(dlx,xv(k,ic)*d**2)                                           2548
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              2549
      o=o+d*(xb(k)/xs(k))                                                  2550
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  2551
16481 continue                                                             2552
16482 continue                                                             2552
      d=0.0                                                                2552
      if(intr.ne.0) d=svr/xm(0)                                            2553
      if(d .eq. 0.0)goto 16531                                             2553
      b(0,ic)=b(0,ic)+d                                                    2553
      dlx=max(dlx,xm(0)*d**2)                                              2554
      r=r-d*v                                                              2554
      svr=svr-d*xm(0)                                                      2555
16531 continue                                                             2556
      if(dlx.lt.shr)goto 16472                                             2556
      if(nlp .le. maxit)goto 16551                                         2556
      jerr=-ilm                                                            2556
      return                                                               2556
16551 continue                                                             2557
      goto 16471                                                           2558
16472 continue                                                             2558
      goto 16341                                                           2559
16342 continue                                                             2559
      if(jxx.gt.0)goto 16252                                               2560
      if(xm(0)*(b(0,ic)-bs(0,ic))**2.gt.shr) ixx=1                         2561
      if(ixx .ne. 0)goto 16571                                             2562
16580 do 16581 j=1,nin                                                     2562
      k=m(j)                                                               2563
      if(xv(k,ic)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 16601                2563
      ixx=1                                                                2563
      goto 16582                                                           2563
16601 continue                                                             2564
16581 continue                                                             2565
16582 continue                                                             2565
16571 continue                                                             2566
      sc=b(0,ic)+g(:,ic)                                                   2566
      b0=0.0                                                               2567
16610 do 16611 j=1,nin                                                     2567
      l=m(j)                                                               2567
      jb=ix(l)                                                             2567
      je=ix(l+1)-1                                                         2568
      sc(jx(jb:je))=sc(jx(jb:je))+b(l,ic)*x(jb:je)/xs(l)                   2569
      b0=b0-b(l,ic)*xb(l)/xs(l)                                            2570
16611 continue                                                             2571
16612 continue                                                             2571
      sc=min(max(exmn,sc+b0),exmx)                                         2572
      sxp=sxp-q(:,ic)                                                      2573
      q(:,ic)=min(max(emin*sxp,exp(sc)),emax*sxp)                          2574
      sxp=sxp+q(:,ic)                                                      2575
16251 continue                                                             2576
16252 continue                                                             2576
      s=-sum(b(0,:))/nc                                                    2576
      b(0,:)=b(0,:)+s                                                      2576
      sc=s                                                                 2576
      b0=0.0                                                               2577
16620 do 16621 j=1,nin                                                     2577
      l=m(j)                                                               2578
      if(vp(l) .gt. 0.0)goto 16641                                         2578
      s=sum(b(l,:))/nc                                                     2578
      goto 16651                                                           2579
16641 continue                                                             2579
      s=elc(parm,nc,cl(:,l),b(l,:),is)                                     2579
16651 continue                                                             2580
16631 continue                                                             2580
      b(l,:)=b(l,:)-s                                                      2581
      jb=ix(l)                                                             2581
      je=ix(l+1)-1                                                         2582
      sc(jx(jb:je))=sc(jx(jb:je))-s*x(jb:je)/xs(l)                         2583
      b0=b0+s*xb(l)/xs(l)                                                  2584
16621 continue                                                             2585
16622 continue                                                             2585
      sc=sc+b0                                                             2585
      sc=exp(sc)                                                           2585
      sxp=sxp*sc                                                           2585
16660 do 16661 ic=1,nc                                                     2585
      q(:,ic)=q(:,ic)*sc                                                   2585
16661 continue                                                             2586
16662 continue                                                             2586
      if(jxx.gt.0)goto 16242                                               2586
      if(ig.eq.0)goto 16242                                                2587
      if(ixx .ne. 0)goto 16681                                             2588
16690 do 16691 j=1,ni                                                      2588
      if(iy(j).eq.1)goto 16691                                             2588
      if(ju(j).eq.0)goto 16691                                             2588
      ga(j)=0.0                                                            2588
16691 continue                                                             2589
16692 continue                                                             2589
16700 do 16701 ic=1,nc                                                     2589
      v=q(:,ic)/sxp                                                        2589
      r=w*(y(:,ic)-v)                                                      2589
      v=w*v*(1.0-v)                                                        2590
16710 do 16711 j=1,ni                                                      2590
      if(iy(j).eq.1)goto 16711                                             2590
      if(ju(j).eq.0)goto 16711                                             2591
      jb=ix(j)                                                             2591
      je=ix(j+1)-1                                                         2591
      jn=ix(j+1)-ix(j)                                                     2592
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2593
      gj=dot_product(sc(1:jn),x(jb:je))                                    2594
      ga(j)=max(ga(j),abs(gj-svr*xb(j))/xs(j))                             2595
16711 continue                                                             2596
16712 continue                                                             2596
16701 continue                                                             2597
16702 continue                                                             2597
16720 do 16721 k=1,ni                                                      2597
      if(iy(k).eq.1)goto 16721                                             2597
      if(ju(k).eq.0)goto 16721                                             2598
      if(ga(k) .le. al1*vp(k))goto 16741                                   2598
      iy(k)=1                                                              2598
      ixx=1                                                                2598
16741 continue                                                             2599
16721 continue                                                             2600
16722 continue                                                             2600
      if(ixx.eq.1) go to 11020                                             2601
      goto 16242                                                           2602
16681 continue                                                             2603
      goto 16241                                                           2604
16242 continue                                                             2604
      if(jxx .le. 0)goto 16761                                             2604
      jerr=-10000-ilm                                                      2604
      goto 16162                                                           2604
16761 continue                                                             2604
      devi=0.0                                                             2605
16770 do 16771 ic=1,nc                                                     2606
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          2606
      a0(ic,ilm)=b(0,ic)                                                   2607
16780 do 16781 i=1,no                                                      2607
      if(y(i,ic).le.0.0)goto 16781                                         2608
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           2609
16781 continue                                                             2610
16782 continue                                                             2610
16771 continue                                                             2611
16772 continue                                                             2611
      kin(ilm)=nin                                                         2611
      alm(ilm)=al                                                          2611
      lmu=ilm                                                              2612
      dev(ilm)=(dev1-devi)/dev0                                            2612
      if(ig.eq.0)goto 16162                                                2613
      if(ilm.lt.mnl)goto 16161                                             2613
      if(flmin.ge.1.0)goto 16161                                           2614
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 16162             2615
      if(dev(ilm).gt.devmax)goto 16162                                     2615
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 16162                             2616
16161 continue                                                             2617
16162 continue                                                             2617
      g=log(q)                                                             2617
16790 do 16791 i=1,no                                                      2617
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         2617
16791 continue                                                             2618
16792 continue                                                             2618
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,xm,sc,ga,iy)                      2619
      return                                                               2620
      end                                                                  2621
      subroutine lcmodval(nc,nx,a0,ca,ia,nin,x,ix,jx,n,f)                  2622
      implicit double precision(a-h,o-z)                                   2623
      double precision a0(nc),ca(nx,nc),x(*),f(nc,n)                       2623
      integer ia(*),ix(*),jx(*)                                            2624
16800 do 16801 ic=1,nc                                                     2624
      f(ic,:)=a0(ic)                                                       2624
16801 continue                                                             2625
16802 continue                                                             2625
16810 do 16811 j=1,nin                                                     2625
      k=ia(j)                                                              2625
      kb=ix(k)                                                             2625
      ke=ix(k+1)-1                                                         2626
16820 do 16821 ic=1,nc                                                     2626
      f(ic,jx(kb:ke))=f(ic,jx(kb:ke))+ca(j,ic)*x(kb:ke)                    2626
16821 continue                                                             2627
16822 continue                                                             2627
16811 continue                                                             2628
16812 continue                                                             2628
      return                                                               2629
      end                                                                  2630
      subroutine coxnet (parm,no,ni,x,y,d,g,w,jd,vp,cl,ne,nx,nlam,flmin,   2632 
     *ulam,thr,  maxit,isd,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   2633
      double precision x(no,ni),y(no),d(no),g(no),w(no),vp(ni),ulam(nlam   2634 
     *)
      double precision ca(nx,nlam),dev(nlam),alm(nlam),cl(2,ni)            2635
      integer jd(*),ia(nx),nin(nlam)                                       2636
      double precision, dimension (:), allocatable :: xs,ww,vq                  
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 16841                                    2640
      jerr=10000                                                           2640
      return                                                               2640
16841 continue                                                             2641
      allocate(ww(1:no),stat=jerr)                                         2642
      if(jerr.ne.0) return                                                 2643
      allocate(ju(1:ni),stat=jerr)                                         2644
      if(jerr.ne.0) return                                                 2645
      allocate(vq(1:ni),stat=jerr)                                         2646
      if(jerr.ne.0) return                                                 2647
      if(isd .le. 0)goto 16861                                             2647
      allocate(xs(1:ni),stat=jerr)                                         2647
      if(jerr.ne.0) return                                                 2647
16861 continue                                                             2649
      call chkvars(no,ni,x,ju)                                             2650
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2651
      if(maxval(ju) .gt. 0)goto 16881                                      2651
      jerr=7777                                                            2651
      return                                                               2651
16881 continue                                                             2652
      vq=max(0d0,vp)                                                       2652
      vq=vq*ni/sum(vq)                                                     2653
      ww=max(0d0,w)                                                        2653
      sw=sum(ww)                                                           2654
      if(sw .gt. 0.0)goto 16901                                            2654
      jerr=9999                                                            2654
      return                                                               2654
16901 continue                                                             2654
      ww=ww/sw                                                             2655
      call cstandard(no,ni,x,ww,ju,isd,xs)                                 2656
      if(isd .le. 0)goto 16921                                             2656
16930 do 16931 j=1,ni                                                      2656
      cl(:,j)=cl(:,j)*xs(j)                                                2656
16931 continue                                                             2656
16932 continue                                                             2656
16921 continue                                                             2657
      call coxnet1(parm,no,ni,x,y,d,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam,   2659 
     *thr,  isd,maxit,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 2659
      dev0=2.0*sw*dev0                                                     2660
      if(isd .le. 0)goto 16951                                             2660
16960 do 16961 k=1,lmu                                                     2660
      nk=nin(k)                                                            2660
      ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                                   2660
16961 continue                                                             2660
16962 continue                                                             2660
16951 continue                                                             2661
      deallocate(ww,ju,vq)                                                 2661
      if(isd.gt.0) deallocate(xs)                                          2662
      return                                                               2663
      end                                                                  2664
      subroutine cstandard (no,ni,x,w,ju,isd,xs)                           2665
      implicit double precision(a-h,o-z)                                   2666
      double precision x(no,ni),w(no),xs(ni)                               2666
      integer ju(ni)                                                       2667
16970 do 16971 j=1,ni                                                      2667
      if(ju(j).eq.0)goto 16971                                             2668
      xm=dot_product(w,x(:,j))                                             2668
      x(:,j)=x(:,j)-xm                                                     2669
      if(isd .le. 0)goto 16991                                             2669
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                                 2669
      x(:,j)=x(:,j)/xs(j)                                                  2669
16991 continue                                                             2670
16971 continue                                                             2671
16972 continue                                                             2671
      return                                                               2672
      end                                                                  2673
      subroutine coxnet1(parm,no,ni,x,y,d,g,q,ju,vp,cl,ne,nx,nlam,flmin,   2675 
     *ulam,cthri,  isd,maxit,lmu,ao,m,kin,dev0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   2676
      double precision x(no,ni),y(no),q(no),d(no),g(no),vp(ni),ulam(nlam   2677 
     *)
      double precision ao(nx,nlam),dev(nlam),alm(nlam),cl(2,ni)            2678
      integer ju(ni),m(nx),kin(nlam)                                       2679
      double precision, dimension (:), allocatable :: w,dk,v,xs,wr              
      double precision, dimension (:), allocatable :: a,as,f,dq                 
      double precision, dimension (:), allocatable :: e,uu,ga                   
      integer, dimension (:), allocatable :: jp,kp,mm,ixx                       
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               2686
      sml=sml*100.0                                                        2686
      devmax=devmax*0.99/0.999                                             2687
      allocate(e(1:no),stat=jerr)                                          2688
      if(jerr.ne.0)go to 12320                                             2689
      allocate(uu(1:no),stat=jerr)                                         2690
      if(jerr.ne.0)go to 12320                                             2691
      allocate(f(1:no),stat=jerr)                                          2692
      if(jerr.ne.0)go to 12320                                             2693
      allocate(w(1:no),stat=jerr)                                          2694
      if(jerr.ne.0)go to 12320                                             2695
      allocate(v(1:ni),stat=jerr)                                          2696
      if(jerr.ne.0)go to 12320                                             2697
      allocate(a(1:ni),stat=jerr)                                          2698
      if(jerr.ne.0)go to 12320                                             2699
      allocate(as(1:ni),stat=jerr)                                         2700
      if(jerr.ne.0)go to 12320                                             2701
      allocate(xs(1:ni),stat=jerr)                                         2702
      if(jerr.ne.0)go to 12320                                             2703
      allocate(ga(1:ni),stat=jerr)                                         2704
      if(jerr.ne.0)go to 12320                                             2705
      allocate(ixx(1:ni),stat=jerr)                                        2706
      if(jerr.ne.0)go to 12320                                             2707
      allocate(jp(1:no),stat=jerr)                                         2708
      if(jerr.ne.0)go to 12320                                             2709
      allocate(kp(1:no),stat=jerr)                                         2710
      if(jerr.ne.0)go to 12320                                             2711
      allocate(dk(1:no),stat=jerr)                                         2712
      if(jerr.ne.0)go to 12320                                             2713
      allocate(wr(1:no),stat=jerr)                                         2714
      if(jerr.ne.0)go to 12320                                             2715
      allocate(dq(1:no),stat=jerr)                                         2716
      if(jerr.ne.0)go to 12320                                             2717
      allocate(mm(1:ni),stat=jerr)                                         2718
      if(jerr.ne.0)go to 12320                                             2719
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                               2720
      if(jerr.ne.0) go to 12320                                            2720
      alpha=parm                                                           2721
      oma=1.0-alpha                                                        2721
      nlm=0                                                                2721
      ixx=0                                                                2721
      al=0.0                                                               2722
      dq=d*q                                                               2722
      call died(no,nk,dq,kp,jp,dk)                                         2723
      a=0.0                                                                2723
      f(1)=0.0                                                             2723
      fmax=log(huge(f(1))*0.1)                                             2724
      if(nonzero(no,g) .eq. 0)goto 17011                                   2724
      f=g-dot_product(q,g)                                                 2725
      e=q*exp(sign(min(abs(f),fmax),f))                                    2726
      goto 17021                                                           2727
17011 continue                                                             2727
      f=0.0                                                                2727
      e=q                                                                  2727
17021 continue                                                             2728
17001 continue                                                             2728
      r0=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                                 2729
      rr=-(dot_product(dk(1:nk),log(dk(1:nk)))+r0)                         2729
      dev0=rr                                                              2730
17030 do 17031 i=1,no                                                      2730
      if((y(i) .ge. t0) .and. (q(i) .gt. 0.0))goto 17051                   2730
      w(i)=0.0                                                             2730
      wr(i)=w(i)                                                           2730
17051 continue                                                             2730
17031 continue                                                             2731
17032 continue                                                             2731
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                         2732
      if(jerr.ne.0) go to 12320                                            2734
      alf=1.0                                                              2736
      if(flmin .ge. 1.0)goto 17071                                         2736
      eqs=max(eps,flmin)                                                   2736
      alf=eqs**(1.0/(nlam-1))                                              2736
17071 continue                                                             2737
      m=0                                                                  2737
      mm=0                                                                 2737
      nlp=0                                                                2737
      nin=nlp                                                              2737
      mnl=min(mnlam,nlam)                                                  2737
      as=0.0                                                               2737
      cthr=cthri*dev0                                                      2738
17080 do 17081 j=1,ni                                                      2738
      if(ju(j).eq.0)goto 17081                                             2738
      ga(j)=abs(dot_product(wr,x(:,j)))                                    2738
17081 continue                                                             2739
17082 continue                                                             2739
17090 do 17091 ilm=1,nlam                                                  2739
      al0=al                                                               2740
      if(flmin .lt. 1.0)goto 17111                                         2740
      al=ulam(ilm)                                                         2740
      goto 17101                                                           2741
17111 if(ilm .le. 2)goto 17121                                             2741
      al=al*alf                                                            2741
      goto 17101                                                           2742
17121 if(ilm .ne. 1)goto 17131                                             2742
      al=big                                                               2742
      goto 17141                                                           2743
17131 continue                                                             2743
      al0=0.0                                                              2744
17150 do 17151 j=1,ni                                                      2744
      if(ju(j).eq.0)goto 17151                                             2744
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2744
17151 continue                                                             2745
17152 continue                                                             2745
      al0=al0/max(parm,1.0d-3)                                             2745
      al=alf*al0                                                           2746
17141 continue                                                             2747
17101 continue                                                             2747
      sa=alpha*al                                                          2747
      omal=oma*al                                                          2747
      tlam=alpha*(2.0*al-al0)                                              2748
17160 do 17161 k=1,ni                                                      2748
      if(ixx(k).eq.1)goto 17161                                            2748
      if(ju(k).eq.0)goto 17161                                             2749
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     2750
17161 continue                                                             2751
17162 continue                                                             2751
11020 continue                                                             2752
17170 continue                                                             2752
17171 continue                                                             2752
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2753
      call vars(no,ni,x,w,ixx,v)                                           2754
17180 continue                                                             2754
17181 continue                                                             2754
      nlp=nlp+1                                                            2754
      dli=0.0                                                              2755
17190 do 17191 j=1,ni                                                      2755
      if(ixx(j).eq.0)goto 17191                                            2756
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                   2757
      if(abs(u) .gt. vp(j)*sa)goto 17211                                   2757
      at=0.0                                                               2757
      goto 17221                                                           2758
17211 continue                                                             2758
      at=max(cl(1,j),min(cl(2,j),sign(abs(u)-vp(j)*sa,u)/  (v(j)+vp(j)*o   2760 
     *mal)))
17221 continue                                                             2761
17201 continue                                                             2761
      if(at .eq. a(j))goto 17241                                           2761
      del=at-a(j)                                                          2761
      a(j)=at                                                              2761
      dli=max(dli,v(j)*del**2)                                             2762
      wr=wr-del*w*x(:,j)                                                   2762
      f=f+del*x(:,j)                                                       2763
      if(mm(j) .ne. 0)goto 17261                                           2763
      nin=nin+1                                                            2763
      if(nin.gt.nx)goto 17192                                              2764
      mm(j)=nin                                                            2764
      m(nin)=j                                                             2765
17261 continue                                                             2766
17241 continue                                                             2767
17191 continue                                                             2768
17192 continue                                                             2768
      if(nin.gt.nx)goto 17182                                              2768
      if(dli.lt.cthr)goto 17182                                            2769
      if(nlp .le. maxit)goto 17281                                         2769
      jerr=-ilm                                                            2769
      return                                                               2769
17281 continue                                                             2770
17290 continue                                                             2770
17291 continue                                                             2770
      nlp=nlp+1                                                            2770
      dli=0.0                                                              2771
17300 do 17301 l=1,nin                                                     2771
      j=m(l)                                                               2772
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                   2773
      if(abs(u) .gt. vp(j)*sa)goto 17321                                   2773
      at=0.0                                                               2773
      goto 17331                                                           2774
17321 continue                                                             2774
      at=max(cl(1,j),min(cl(2,j),sign(abs(u)-vp(j)*sa,u)/  (v(j)+vp(j)*o   2776 
     *mal)))
17331 continue                                                             2777
17311 continue                                                             2777
      if(at .eq. a(j))goto 17351                                           2777
      del=at-a(j)                                                          2777
      a(j)=at                                                              2777
      dli=max(dli,v(j)*del**2)                                             2778
      wr=wr-del*w*x(:,j)                                                   2778
      f=f+del*x(:,j)                                                       2779
17351 continue                                                             2780
17301 continue                                                             2781
17302 continue                                                             2781
      if(dli.lt.cthr)goto 17292                                            2781
      if(nlp .le. maxit)goto 17371                                         2781
      jerr=-ilm                                                            2781
      return                                                               2781
17371 continue                                                             2782
      goto 17291                                                           2783
17292 continue                                                             2783
      goto 17181                                                           2784
17182 continue                                                             2784
      if(nin.gt.nx)goto 17172                                              2785
      e=q*exp(sign(min(abs(f),fmax),f))                                    2786
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                         2787
      if(jerr .eq. 0)goto 17391                                            2787
      jerr=jerr-ilm                                                        2787
      go to 12320                                                          2787
17391 continue                                                             2788
      ix=0                                                                 2789
17400 do 17401 j=1,nin                                                     2789
      k=m(j)                                                               2790
      if(v(k)*(a(k)-as(k))**2.lt.cthr)goto 17401                           2790
      ix=1                                                                 2790
      goto 17402                                                           2790
17401 continue                                                             2791
17402 continue                                                             2791
      if(ix .ne. 0)goto 17421                                              2792
17430 do 17431 k=1,ni                                                      2792
      if(ixx(k).eq.1)goto 17431                                            2792
      if(ju(k).eq.0)goto 17431                                             2793
      ga(k)=abs(dot_product(wr,x(:,k)))                                    2794
      if(ga(k) .le. sa*vp(k))goto 17451                                    2794
      ixx(k)=1                                                             2794
      ix=1                                                                 2794
17451 continue                                                             2795
17431 continue                                                             2796
17432 continue                                                             2796
      if(ix.eq.1) go to 11020                                              2797
      goto 17172                                                           2798
17421 continue                                                             2799
      goto 17171                                                           2800
17172 continue                                                             2800
      if(nin .le. nx)goto 17471                                            2800
      jerr=-10000-ilm                                                      2800
      goto 17092                                                           2800
17471 continue                                                             2801
      if(nin.gt.0) ao(1:nin,ilm)=a(m(1:nin))                               2801
      kin(ilm)=nin                                                         2802
      alm(ilm)=al                                                          2802
      lmu=ilm                                                              2803
      dev(ilm)=(risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)-r0)/rr                   2804
      if(ilm.lt.mnl)goto 17091                                             2804
      if(flmin.ge.1.0)goto 17091                                           2805
      me=0                                                                 2805
17480 do 17481 j=1,nin                                                     2805
      if(ao(j,ilm).ne.0.0) me=me+1                                         2805
17481 continue                                                             2805
17482 continue                                                             2805
      if(me.gt.ne)goto 17092                                               2806
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 17092              2807
      if(dev(ilm).gt.devmax)goto 17092                                     2808
17091 continue                                                             2809
17092 continue                                                             2809
      g=f                                                                  2810
12320 continue                                                             2810
      deallocate(e,uu,w,dk,v,xs,f,wr,a,as,jp,kp,dq,mm,ga,ixx)              2811
      return                                                               2812
      end                                                                  2813
      subroutine cxmodval(ca,ia,nin,n,x,f)                                 2814
      implicit double precision(a-h,o-z)                                   2815
      double precision ca(nin),x(n,*),f(n)                                 2815
      integer ia(nin)                                                      2816
      f=0.0                                                                2816
      if(nin.le.0) return                                                  2817
17490 do 17491 i=1,n                                                       2817
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                      2817
17491 continue                                                             2818
17492 continue                                                             2818
      return                                                               2819
      end                                                                  2820
      subroutine groups(no,y,d,q,nk,kp,jp,t0,jerr)                         2821
      implicit double precision(a-h,o-z)                                   2822
      double precision y(no),d(no),q(no)                                   2822
      integer jp(no),kp(*)                                                 2823
17500 do 17501 j=1,no                                                      2823
      jp(j)=j                                                              2823
17501 continue                                                             2823
17502 continue                                                             2823
      call psort7(y,jp,1,no)                                               2824
      nj=0                                                                 2824
17510 do 17511 j=1,no                                                      2824
      if(q(jp(j)).le.0.0)goto 17511                                        2824
      nj=nj+1                                                              2824
      jp(nj)=jp(j)                                                         2824
17511 continue                                                             2825
17512 continue                                                             2825
      if(nj .ne. 0)goto 17531                                              2825
      jerr=20000                                                           2825
      return                                                               2825
17531 continue                                                             2826
      j=1                                                                  2826
17540 continue                                                             2826
17541 if(d(jp(j)).gt.0.0)goto 17542                                        2826
      j=j+1                                                                2826
      if(j.gt.nj)goto 17542                                                2826
      goto 17541                                                           2827
17542 continue                                                             2827
      if(j .lt. nj-1)goto 17561                                            2827
      jerr=30000                                                           2827
      return                                                               2827
17561 continue                                                             2828
      t0=y(jp(j))                                                          2828
      j0=j-1                                                               2829
      if(j0 .le. 0)goto 17581                                              2830
17590 continue                                                             2830
17591 if(y(jp(j0)).lt.t0)goto 17592                                        2830
      j0=j0-1                                                              2830
      if(j0.eq.0)goto 17592                                                2830
      goto 17591                                                           2831
17592 continue                                                             2831
      if(j0 .le. 0)goto 17611                                              2831
      nj=nj-j0                                                             2831
17620 do 17621 j=1,nj                                                      2831
      jp(j)=jp(j+j0)                                                       2831
17621 continue                                                             2831
17622 continue                                                             2831
17611 continue                                                             2832
17581 continue                                                             2833
      jerr=0                                                               2833
      nk=0                                                                 2833
      yk=t0                                                                2833
      j=2                                                                  2834
17630 continue                                                             2834
17631 continue                                                             2834
17640 continue                                                             2835
17641 if(d(jp(j)).gt.0.0.and.y(jp(j)).gt.yk)goto 17642                     2835
      j=j+1                                                                2835
      if(j.gt.nj)goto 17642                                                2835
      goto 17641                                                           2836
17642 continue                                                             2836
      nk=nk+1                                                              2836
      kp(nk)=j-1                                                           2836
      if(j.gt.nj)goto 17632                                                2837
      if(j .ne. nj)goto 17661                                              2837
      nk=nk+1                                                              2837
      kp(nk)=nj                                                            2837
      goto 17632                                                           2837
17661 continue                                                             2838
      yk=y(jp(j))                                                          2838
      j=j+1                                                                2839
      goto 17631                                                           2840
17632 continue                                                             2840
      return                                                               2841
      end                                                                  2842
      subroutine outer(no,nk,d,dk,kp,jp,e,wr,w,jerr,u)                     2843
      implicit double precision(a-h,o-z)                                   2844
      double precision d(no),dk(nk),wr(no),w(no)                           2845
      double precision e(no),u(no),b,c                                     2845
      integer kp(nk),jp(no)                                                2846
      call usk(no,nk,kp,jp,e,u)                                            2847
      b=dk(1)/u(1)                                                         2847
      c=dk(1)/u(1)**2                                                      2847
      jerr=0                                                               2848
17670 do 17671 j=1,kp(1)                                                   2848
      i=jp(j)                                                              2849
      w(i)=e(i)*(b-e(i)*c)                                                 2849
      if(w(i) .gt. 0.0)goto 17691                                          2849
      jerr=-30000                                                          2849
      return                                                               2849
17691 continue                                                             2850
      wr(i)=d(i)-e(i)*b                                                    2851
17671 continue                                                             2852
17672 continue                                                             2852
17700 do 17701 k=2,nk                                                      2852
      j1=kp(k-1)+1                                                         2852
      j2=kp(k)                                                             2853
      b=b+dk(k)/u(k)                                                       2853
      c=c+dk(k)/u(k)**2                                                    2854
17710 do 17711 j=j1,j2                                                     2854
      i=jp(j)                                                              2855
      w(i)=e(i)*(b-e(i)*c)                                                 2855
      if(w(i) .gt. 0.0)goto 17731                                          2855
      jerr=-30000                                                          2855
      return                                                               2855
17731 continue                                                             2856
      wr(i)=d(i)-e(i)*b                                                    2857
17711 continue                                                             2858
17712 continue                                                             2858
17701 continue                                                             2859
17702 continue                                                             2859
      return                                                               2860
      end                                                                  2861
      subroutine vars(no,ni,x,w,ixx,v)                                     2862
      implicit double precision(a-h,o-z)                                   2863
      double precision x(no,ni),w(no),v(ni)                                2863
      integer ixx(ni)                                                      2864
17740 do 17741 j=1,ni                                                      2864
      if(ixx(j).gt.0) v(j)=dot_product(w,x(:,j)**2)                        2864
17741 continue                                                             2865
17742 continue                                                             2865
      return                                                               2866
      end                                                                  2867
      subroutine died(no,nk,d,kp,jp,dk)                                    2868
      implicit double precision(a-h,o-z)                                   2869
      double precision d(no),dk(nk)                                        2869
      integer kp(nk),jp(no)                                                2870
      dk(1)=sum(d(jp(1:kp(1))))                                            2871
17750 do 17751 k=2,nk                                                      2871
      dk(k)=sum(d(jp((kp(k-1)+1):kp(k))))                                  2871
17751 continue                                                             2872
17752 continue                                                             2872
      return                                                               2873
      end                                                                  2874
      subroutine usk(no,nk,kp,jp,e,u)                                      2875
      implicit double precision(a-h,o-z)                                   2876
      double precision e(no),u(nk),h                                       2876
      integer kp(nk),jp(no)                                                2877
      h=0.0                                                                2878
17760 do 17761 k=nk,1,-1                                                   2878
      j2=kp(k)                                                             2879
      j1=1                                                                 2879
      if(k.gt.1) j1=kp(k-1)+1                                              2880
17770 do 17771 j=j2,j1,-1                                                  2880
      h=h+e(jp(j))                                                         2880
17771 continue                                                             2881
17772 continue                                                             2881
      u(k)=h                                                               2882
17761 continue                                                             2883
17762 continue                                                             2883
      return                                                               2884
      end                                                                  2885
      function risk(no,ni,nk,d,dk,f,e,kp,jp,u)                             2886
      implicit double precision(a-h,o-z)                                   2887
      double precision d(no),dk(nk),f(no)                                  2888
      integer kp(nk),jp(no)                                                2888
      double precision e(no),u(nk),s                                       2889
      call usk(no,nk,kp,jp,e,u)                                            2889
      u=log(u)                                                             2890
      risk=dot_product(d,f)-dot_product(dk,u)                              2891
      return                                                               2892
      end                                                                  2893
      subroutine loglike(no,ni,x,y,d,g,w,nlam,a,flog,jerr)                 2894
      implicit double precision(a-h,o-z)                                   2895
      double precision x(no,ni),y(no),d(no),g(no),w(no),a(ni,nlam),flog(   2896 
     *nlam)
      double precision, dimension (:), allocatable :: dk,f,xm,dq,q              
      double precision, dimension (:), allocatable :: e,uu                      
      integer, dimension (:), allocatable :: jp,kp                              
      allocate(e(1:no),stat=jerr)                                          2902
      if(jerr.ne.0) go to 12320                                            2903
      allocate(q(1:no),stat=jerr)                                          2904
      if(jerr.ne.0) go to 12320                                            2905
      allocate(uu(1:no),stat=jerr)                                         2906
      if(jerr.ne.0) go to 12320                                            2907
      allocate(f(1:no),stat=jerr)                                          2908
      if(jerr.ne.0) go to 12320                                            2909
      allocate(dk(1:no),stat=jerr)                                         2910
      if(jerr.ne.0) go to 12320                                            2911
      allocate(jp(1:no),stat=jerr)                                         2912
      if(jerr.ne.0) go to 12320                                            2913
      allocate(kp(1:no),stat=jerr)                                         2914
      if(jerr.ne.0) go to 12320                                            2915
      allocate(dq(1:no),stat=jerr)                                         2916
      if(jerr.ne.0) go to 12320                                            2917
      allocate(xm(1:ni),stat=jerr)                                         2918
      if(jerr.ne.0) go to 12320                                            2919
      q=max(0d0,w)                                                         2919
      sw=sum(q)                                                            2920
      if(sw .gt. 0.0)goto 17791                                            2920
      jerr=9999                                                            2920
      go to 12320                                                          2920
17791 continue                                                             2921
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                               2922
      if(jerr.ne.0) go to 12320                                            2922
      fmax=log(huge(e(1))*0.1)                                             2923
      dq=d*q                                                               2923
      call died(no,nk,dq,kp,jp,dk)                                         2923
      gm=dot_product(q,g)/sw                                               2924
17800 do 17801 j=1,ni                                                      2924
      xm(j)=dot_product(q,x(:,j))/sw                                       2924
17801 continue                                                             2925
17802 continue                                                             2925
17810 do 17811 lam=1,nlam                                                  2926
17820 do 17821 i=1,no                                                      2926
      f(i)=g(i)-gm+dot_product(a(:,lam),(x(i,:)-xm))                       2927
      e(i)=q(i)*exp(sign(min(abs(f(i)),fmax),f(i)))                        2928
17821 continue                                                             2929
17822 continue                                                             2929
      flog(lam)=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                          2930
17811 continue                                                             2931
17812 continue                                                             2931
12320 continue                                                             2931
      deallocate(e,uu,dk,f,jp,kp,dq)                                       2932
      return                                                               2933
      end                                                                  2934
      subroutine fishnet (parm,no,ni,x,y,g,w,jd,vp,cl,ne,nx,nlam,flmin,u   2936 
     *lam,thr,  isd,intr,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   2937
      double precision x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam)        2938
      double precision ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam),cl(2,ni)   2939
      integer jd(*),ia(nx),nin(nlam)                                       2940
      double precision, dimension (:), allocatable :: xm,xs,ww,vq               
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 17841                                    2944
      jerr=10000                                                           2944
      return                                                               2944
17841 continue                                                             2945
      if(minval(y) .ge. 0.0)goto 17861                                     2945
      jerr=8888                                                            2945
      return                                                               2945
17861 continue                                                             2946
      allocate(ww(1:no),stat=jerr)                                         2947
      if(jerr.ne.0) return                                                 2948
      allocate(ju(1:ni),stat=jerr)                                         2949
      if(jerr.ne.0) return                                                 2950
      allocate(vq(1:ni),stat=jerr)                                         2951
      if(jerr.ne.0) return                                                 2952
      allocate(xm(1:ni),stat=jerr)                                         2953
      if(jerr.ne.0) return                                                 2954
      if(isd .le. 0)goto 17881                                             2954
      allocate(xs(1:ni),stat=jerr)                                         2954
      if(jerr.ne.0) return                                                 2954
17881 continue                                                             2955
      call chkvars(no,ni,x,ju)                                             2956
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2957
      if(maxval(ju) .gt. 0)goto 17901                                      2957
      jerr=7777                                                            2957
      go to 12320                                                          2957
17901 continue                                                             2958
      vq=max(0d0,vp)                                                       2958
      vq=vq*ni/sum(vq)                                                     2959
      ww=max(0d0,w)                                                        2959
      sw=sum(ww)                                                           2959
      if(sw .gt. 0.0)goto 17921                                            2959
      jerr=9999                                                            2959
      go to 12320                                                          2959
17921 continue                                                             2960
      ww=ww/sw                                                             2961
      call lstandard1(no,ni,x,ww,ju,isd,intr,xm,xs)                        2962
      if(isd .le. 0)goto 17941                                             2962
17950 do 17951 j=1,ni                                                      2962
      cl(:,j)=cl(:,j)*xs(j)                                                2962
17951 continue                                                             2962
17952 continue                                                             2962
17941 continue                                                             2963
      call fishnet1(parm,no,ni,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam,t   2965 
     *hr,  isd,intr,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) go to 12320                                            2965
      dev0=2.0*sw*dev0                                                     2966
17960 do 17961 k=1,lmu                                                     2966
      nk=nin(k)                                                            2967
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                      2968
      if(intr .ne. 0)goto 17981                                            2968
      a0(k)=0.0                                                            2968
      goto 17991                                                           2969
17981 continue                                                             2969
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                     2969
17991 continue                                                             2970
17971 continue                                                             2970
17961 continue                                                             2971
17962 continue                                                             2971
12320 continue                                                             2971
      deallocate(ww,ju,vq,xm)                                              2971
      if(isd.gt.0) deallocate(xs)                                          2972
      return                                                               2973
      end                                                                  2974
      subroutine fishnet1(parm,no,ni,x,y,g,q,ju,vp,cl,ne,nx,nlam,flmin,u   2976 
     *lam,shri,  isd,intr,maxit,lmu,a0,ca,m,kin,dev0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   2977
      double precision x(no,ni),y(no),g(no),q(no),vp(ni),ulam(nlam)        2978
      double precision ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam),cl(2,ni)   2979
      integer ju(ni),m(nx),kin(nlam)                                       2980
      double precision, dimension (:), allocatable :: t,w,wr,v,a,f,as,ga        
      integer, dimension (:), allocatable :: mm,ixx                             
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               2984
      sml=sml*10.0                                                         2985
      allocate(a(1:ni),stat=jerr)                                          2986
      if(jerr.ne.0) return                                                 2987
      allocate(as(1:ni),stat=jerr)                                         2988
      if(jerr.ne.0) return                                                 2989
      allocate(t(1:no),stat=jerr)                                          2990
      if(jerr.ne.0) return                                                 2991
      allocate(mm(1:ni),stat=jerr)                                         2992
      if(jerr.ne.0) return                                                 2993
      allocate(ga(1:ni),stat=jerr)                                         2994
      if(jerr.ne.0) return                                                 2995
      allocate(ixx(1:ni),stat=jerr)                                        2996
      if(jerr.ne.0) return                                                 2997
      allocate(wr(1:no),stat=jerr)                                         2998
      if(jerr.ne.0) return                                                 2999
      allocate(v(1:ni),stat=jerr)                                          3000
      if(jerr.ne.0) return                                                 3001
      allocate(w(1:no),stat=jerr)                                          3002
      if(jerr.ne.0) return                                                 3003
      allocate(f(1:no),stat=jerr)                                          3004
      if(jerr.ne.0) return                                                 3005
      bta=parm                                                             3005
      omb=1.0-bta                                                          3006
      t=q*y                                                                3006
      yb=sum(t)                                                            3006
      fmax=log(huge(bta)*0.1)                                              3007
      if(nonzero(no,g) .ne. 0)goto 18011                                   3008
      if(intr .eq. 0)goto 18031                                            3008
      w=q*yb                                                               3008
      az=log(yb)                                                           3008
      f=az                                                                 3008
      dv0=yb*(az-1.0)                                                      3008
      goto 18041                                                           3009
18031 continue                                                             3009
      w=q                                                                  3009
      az=0.0                                                               3009
      f=az                                                                 3009
      dv0=-1.0                                                             3009
18041 continue                                                             3010
18021 continue                                                             3010
      goto 18051                                                           3011
18011 continue                                                             3011
      w=q*exp(sign(min(abs(g),fmax),g))                                    3011
      v0=sum(w)                                                            3012
      if(intr .eq. 0)goto 18071                                            3012
      eaz=yb/v0                                                            3012
      w=eaz*w                                                              3012
      az=log(eaz)                                                          3012
      f=az+g                                                               3013
      dv0=dot_product(t,g)-yb*(1.0-az)                                     3014
      goto 18081                                                           3015
18071 continue                                                             3015
      az=0.0                                                               3015
      f=g                                                                  3015
      dv0=dot_product(t,g)-v0                                              3015
18081 continue                                                             3016
18061 continue                                                             3016
18051 continue                                                             3017
18001 continue                                                             3017
      a=0.0                                                                3017
      as=0.0                                                               3017
      wr=t-w                                                               3017
      v0=1.0                                                               3017
      if(intr.ne.0) v0=yb                                                  3017
      dvr=-yb                                                              3018
18090 do 18091 i=1,no                                                      3018
      if(t(i).gt.0.0) dvr=dvr+t(i)*log(y(i))                               3018
18091 continue                                                             3018
18092 continue                                                             3018
      dvr=dvr-dv0                                                          3018
      dev0=dvr                                                             3020
      alf=1.0                                                              3022
      if(flmin .ge. 1.0)goto 18111                                         3022
      eqs=max(eps,flmin)                                                   3022
      alf=eqs**(1.0/(nlam-1))                                              3022
18111 continue                                                             3023
      m=0                                                                  3023
      mm=0                                                                 3023
      nlp=0                                                                3023
      nin=nlp                                                              3023
      mnl=min(mnlam,nlam)                                                  3023
      shr=shri*dev0                                                        3023
      ixx=0                                                                3023
      al=0.0                                                               3024
18120 do 18121 j=1,ni                                                      3024
      if(ju(j).eq.0)goto 18121                                             3024
      ga(j)=abs(dot_product(wr,x(:,j)))                                    3024
18121 continue                                                             3025
18122 continue                                                             3025
18130 do 18131 ilm=1,nlam                                                  3025
      al0=al                                                               3026
      if(flmin .lt. 1.0)goto 18151                                         3026
      al=ulam(ilm)                                                         3026
      goto 18141                                                           3027
18151 if(ilm .le. 2)goto 18161                                             3027
      al=al*alf                                                            3027
      goto 18141                                                           3028
18161 if(ilm .ne. 1)goto 18171                                             3028
      al=big                                                               3028
      goto 18181                                                           3029
18171 continue                                                             3029
      al0=0.0                                                              3030
18190 do 18191 j=1,ni                                                      3030
      if(ju(j).eq.0)goto 18191                                             3030
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            3030
18191 continue                                                             3031
18192 continue                                                             3031
      al0=al0/max(bta,1.0d-3)                                              3031
      al=alf*al0                                                           3032
18181 continue                                                             3033
18141 continue                                                             3033
      al2=al*omb                                                           3033
      al1=al*bta                                                           3033
      tlam=bta*(2.0*al-al0)                                                3034
18200 do 18201 k=1,ni                                                      3034
      if(ixx(k).eq.1)goto 18201                                            3034
      if(ju(k).eq.0)goto 18201                                             3035
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     3036
18201 continue                                                             3037
18202 continue                                                             3037
11020 continue                                                             3038
18210 continue                                                             3038
18211 continue                                                             3038
      az0=az                                                               3039
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                3040
18220 do 18221 j=1,ni                                                      3040
      if(ixx(j).ne.0) v(j)=dot_product(w,x(:,j)**2)                        3040
18221 continue                                                             3041
18222 continue                                                             3041
18230 continue                                                             3041
18231 continue                                                             3041
      nlp=nlp+1                                                            3041
      dlx=0.0                                                              3042
18240 do 18241 k=1,ni                                                      3042
      if(ixx(k).eq.0)goto 18241                                            3042
      ak=a(k)                                                              3043
      u=dot_product(wr,x(:,k))+v(k)*ak                                     3043
      au=abs(u)-vp(k)*al1                                                  3044
      if(au .gt. 0.0)goto 18261                                            3044
      a(k)=0.0                                                             3044
      goto 18271                                                           3045
18261 continue                                                             3046
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))           3047
18271 continue                                                             3048
18251 continue                                                             3048
      if(a(k).eq.ak)goto 18241                                             3048
      d=a(k)-ak                                                            3048
      dlx=max(dlx,v(k)*d**2)                                               3049
      wr=wr-d*w*x(:,k)                                                     3049
      f=f+d*x(:,k)                                                         3050
      if(mm(k) .ne. 0)goto 18291                                           3050
      nin=nin+1                                                            3050
      if(nin.gt.nx)goto 18242                                              3051
      mm(k)=nin                                                            3051
      m(nin)=k                                                             3052
18291 continue                                                             3053
18241 continue                                                             3054
18242 continue                                                             3054
      if(nin.gt.nx)goto 18232                                              3055
      if(intr .eq. 0)goto 18311                                            3055
      d=sum(wr)/v0                                                         3056
      az=az+d                                                              3056
      dlx=max(dlx,v0*d**2)                                                 3056
      wr=wr-d*w                                                            3056
      f=f+d                                                                3057
18311 continue                                                             3058
      if(dlx.lt.shr)goto 18232                                             3058
      if(nlp .le. maxit)goto 18331                                         3058
      jerr=-ilm                                                            3058
      return                                                               3058
18331 continue                                                             3059
18340 continue                                                             3059
18341 continue                                                             3059
      nlp=nlp+1                                                            3059
      dlx=0.0                                                              3060
18350 do 18351 l=1,nin                                                     3060
      k=m(l)                                                               3060
      ak=a(k)                                                              3061
      u=dot_product(wr,x(:,k))+v(k)*ak                                     3061
      au=abs(u)-vp(k)*al1                                                  3062
      if(au .gt. 0.0)goto 18371                                            3062
      a(k)=0.0                                                             3062
      goto 18381                                                           3063
18371 continue                                                             3064
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))           3065
18381 continue                                                             3066
18361 continue                                                             3066
      if(a(k).eq.ak)goto 18351                                             3066
      d=a(k)-ak                                                            3066
      dlx=max(dlx,v(k)*d**2)                                               3067
      wr=wr-d*w*x(:,k)                                                     3067
      f=f+d*x(:,k)                                                         3069
18351 continue                                                             3069
18352 continue                                                             3069
      if(intr .eq. 0)goto 18401                                            3069
      d=sum(wr)/v0                                                         3069
      az=az+d                                                              3070
      dlx=max(dlx,v0*d**2)                                                 3070
      wr=wr-d*w                                                            3070
      f=f+d                                                                3071
18401 continue                                                             3072
      if(dlx.lt.shr)goto 18342                                             3072
      if(nlp .le. maxit)goto 18421                                         3072
      jerr=-ilm                                                            3072
      return                                                               3072
18421 continue                                                             3073
      goto 18341                                                           3074
18342 continue                                                             3074
      goto 18231                                                           3075
18232 continue                                                             3075
      if(nin.gt.nx)goto 18212                                              3076
      w=q*exp(sign(min(abs(f),fmax),f))                                    3076
      v0=sum(w)                                                            3076
      wr=t-w                                                               3077
      if(v0*(az-az0)**2 .ge. shr)goto 18441                                3077
      ix=0                                                                 3078
18450 do 18451 j=1,nin                                                     3078
      k=m(j)                                                               3079
      if(v(k)*(a(k)-as(k))**2.lt.shr)goto 18451                            3079
      ix=1                                                                 3079
      goto 18452                                                           3080
18451 continue                                                             3081
18452 continue                                                             3081
      if(ix .ne. 0)goto 18471                                              3082
18480 do 18481 k=1,ni                                                      3082
      if(ixx(k).eq.1)goto 18481                                            3082
      if(ju(k).eq.0)goto 18481                                             3083
      ga(k)=abs(dot_product(wr,x(:,k)))                                    3084
      if(ga(k) .le. al1*vp(k))goto 18501                                   3084
      ixx(k)=1                                                             3084
      ix=1                                                                 3084
18501 continue                                                             3085
18481 continue                                                             3086
18482 continue                                                             3086
      if(ix.eq.1) go to 11020                                              3087
      goto 18212                                                           3088
18471 continue                                                             3089
18441 continue                                                             3090
      goto 18211                                                           3091
18212 continue                                                             3091
      if(nin .le. nx)goto 18521                                            3091
      jerr=-10000-ilm                                                      3091
      goto 18132                                                           3091
18521 continue                                                             3092
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                               3092
      kin(ilm)=nin                                                         3093
      a0(ilm)=az                                                           3093
      alm(ilm)=al                                                          3093
      lmu=ilm                                                              3094
      dev(ilm)=(dot_product(t,f)-v0-dv0)/dvr                               3095
      if(ilm.lt.mnl)goto 18131                                             3095
      if(flmin.ge.1.0)goto 18131                                           3096
      me=0                                                                 3096
18530 do 18531 j=1,nin                                                     3096
      if(ca(j,ilm).ne.0.0) me=me+1                                         3096
18531 continue                                                             3096
18532 continue                                                             3096
      if(me.gt.ne)goto 18132                                               3097
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 18132              3098
      if(dev(ilm).gt.devmax)goto 18132                                     3099
18131 continue                                                             3100
18132 continue                                                             3100
      g=f                                                                  3101
12320 continue                                                             3101
      deallocate(t,w,wr,v,a,f,as,mm,ga,ixx)                                3102
      return                                                               3103
      end                                                                  3104
      function nonzero(n,v)                                                3105
      implicit double precision(a-h,o-z)                                   3106
      double precision v(n)                                                3107
      nonzero=0                                                            3107
18540 do 18541 i=1,n                                                       3107
      if(v(i) .eq. 0.0)goto 18561                                          3107
      nonzero=1                                                            3107
      return                                                               3107
18561 continue                                                             3107
18541 continue                                                             3108
18542 continue                                                             3108
      return                                                               3109
      end                                                                  3110
      subroutine solns(ni,nx,lmu,a,ia,nin,b)                               3111
      implicit double precision(a-h,o-z)                                   3112
      double precision a(nx,lmu),b(ni,lmu)                                 3112
      integer ia(nx),nin(lmu)                                              3113
18570 do 18571 lam=1,lmu                                                   3113
      call uncomp(ni,a(:,lam),ia,nin(lam),b(:,lam))                        3113
18571 continue                                                             3114
18572 continue                                                             3114
      return                                                               3115
      end                                                                  3116
      subroutine lsolns(ni,nx,nc,lmu,a,ia,nin,b)                           3117
      implicit double precision(a-h,o-z)                                   3118
      double precision a(nx,nc,lmu),b(ni,nc,lmu)                           3118
      integer ia(nx),nin(lmu)                                              3119
18580 do 18581 lam=1,lmu                                                   3119
      call luncomp(ni,nx,nc,a(1,1,lam),ia,nin(lam),b(1,1,lam))             3119
18581 continue                                                             3120
18582 continue                                                             3120
      return                                                               3121
      end                                                                  3122
      subroutine deviance(no,ni,x,y,g,q,nlam,a0,a,flog,jerr)               3123
      implicit double precision(a-h,o-z)                                   3124
      double precision x(no,ni),y(no),g(no),q(no),a(ni,nlam),a0(nlam),fl   3125 
     *og(nlam)
      double precision, dimension (:), allocatable :: w                         
      if(minval(y) .ge. 0.0)goto 18601                                     3128
      jerr=8888                                                            3128
      return                                                               3128
18601 continue                                                             3129
      allocate(w(1:no),stat=jerr)                                          3129
      if(jerr.ne.0) return                                                 3130
      w=max(0d0,q)                                                         3130
      sw=sum(w)                                                            3130
      if(sw .gt. 0.0)goto 18621                                            3130
      jerr=9999                                                            3130
      go to 12320                                                          3130
18621 continue                                                             3131
      yb=dot_product(w,y)/sw                                               3131
      fmax=log(huge(y(1))*0.1)                                             3132
18630 do 18631 lam=1,nlam                                                  3132
      s=0.0                                                                3133
18640 do 18641 i=1,no                                                      3133
      if(w(i).le.0.0)goto 18641                                            3134
      f=g(i)+a0(lam)+dot_product(a(:,lam),x(i,:))                          3135
      s=s+w(i)*(y(i)*f-exp(sign(min(abs(f),fmax),f)))                      3136
18641 continue                                                             3137
18642 continue                                                             3137
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                3138
18631 continue                                                             3139
18632 continue                                                             3139
12320 continue                                                             3139
      deallocate(w)                                                        3140
      return                                                               3141
      end                                                                  3142
      subroutine spfishnet (parm,no,ni,x,ix,jx,y,g,w,jd,vp,cl,ne,nx,nlam   3144 
     *,flmin,  ulam,thr,isd,intr,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp
     *,jerr)
      implicit double precision(a-h,o-z)                                   3145
      double precision x(*),y(no),g(no),w(no),vp(ni),ulam(nlam),cl(2,ni)   3146
      double precision ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)            3147
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           3148
      double precision, dimension (:), allocatable :: xm,xs,ww,vq               
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 18661                                    3152
      jerr=10000                                                           3152
      return                                                               3152
18661 continue                                                             3153
      if(minval(y) .ge. 0.0)goto 18681                                     3153
      jerr=8888                                                            3153
      return                                                               3153
18681 continue                                                             3154
      allocate(ww(1:no),stat=jerr)                                         3155
      if(jerr.ne.0) return                                                 3156
      allocate(ju(1:ni),stat=jerr)                                         3157
      if(jerr.ne.0) return                                                 3158
      allocate(vq(1:ni),stat=jerr)                                         3159
      if(jerr.ne.0) return                                                 3160
      allocate(xm(1:ni),stat=jerr)                                         3161
      if(jerr.ne.0) return                                                 3162
      allocate(xs(1:ni),stat=jerr)                                         3163
      if(jerr.ne.0) return                                                 3164
      call spchkvars(no,ni,x,ix,ju)                                        3165
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 3166
      if(maxval(ju) .gt. 0)goto 18701                                      3166
      jerr=7777                                                            3166
      go to 12320                                                          3166
18701 continue                                                             3167
      vq=max(0d0,vp)                                                       3167
      vq=vq*ni/sum(vq)                                                     3168
      ww=max(0d0,w)                                                        3168
      sw=sum(ww)                                                           3168
      if(sw .gt. 0.0)goto 18721                                            3168
      jerr=9999                                                            3168
      go to 12320                                                          3168
18721 continue                                                             3169
      ww=ww/sw                                                             3170
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs)                3171
      if(isd .le. 0)goto 18741                                             3171
18750 do 18751 j=1,ni                                                      3171
      cl(:,j)=cl(:,j)*xs(j)                                                3171
18751 continue                                                             3171
18752 continue                                                             3171
18741 continue                                                             3172
      call spfishnet1(parm,no,ni,x,ix,jx,y,g,ww,ju,vq,cl,ne,nx,nlam,flmi   3174 
     *n,ulam,thr,  isd,intr,maxit,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm,nl
     *p,jerr)
      if(jerr.gt.0) go to 12320                                            3174
      dev0=2.0*sw*dev0                                                     3175
18760 do 18761 k=1,lmu                                                     3175
      nk=nin(k)                                                            3176
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                      3177
      if(intr .ne. 0)goto 18781                                            3177
      a0(k)=0.0                                                            3177
      goto 18791                                                           3178
18781 continue                                                             3178
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                     3178
18791 continue                                                             3179
18771 continue                                                             3179
18761 continue                                                             3180
18762 continue                                                             3180
12320 continue                                                             3180
      deallocate(ww,ju,vq,xm,xs)                                           3181
      return                                                               3182
      end                                                                  3183
      subroutine spfishnet1(parm,no,ni,x,ix,jx,y,g,q,ju,vp,cl,ne,nx,nlam   3185 
     *,flmin,ulam,  shri,isd,intr,maxit,xb,xs,lmu,a0,ca,m,kin,dev0,dev,a
     *lm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   3186
      double precision x(*),y(no),g(no),q(no),vp(ni),ulam(nlam),xb(ni),x   3187 
     *s(ni)
      double precision ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam),cl(2,ni)   3188
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           3189
      double precision, dimension (:), allocatable :: qy,t,w,wr,v               
      double precision, dimension (:), allocatable :: a,as,xm,ga                
      integer, dimension (:), allocatable :: mm,ixx                             
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               3194
      sml=sml*10.0                                                         3195
      allocate(a(1:ni),stat=jerr)                                          3196
      if(jerr.ne.0) return                                                 3197
      allocate(as(1:ni),stat=jerr)                                         3198
      if(jerr.ne.0) return                                                 3199
      allocate(t(1:no),stat=jerr)                                          3200
      if(jerr.ne.0) return                                                 3201
      allocate(mm(1:ni),stat=jerr)                                         3202
      if(jerr.ne.0) return                                                 3203
      allocate(ga(1:ni),stat=jerr)                                         3204
      if(jerr.ne.0) return                                                 3205
      allocate(ixx(1:ni),stat=jerr)                                        3206
      if(jerr.ne.0) return                                                 3207
      allocate(wr(1:no),stat=jerr)                                         3208
      if(jerr.ne.0) return                                                 3209
      allocate(v(1:ni),stat=jerr)                                          3210
      if(jerr.ne.0) return                                                 3211
      allocate(xm(1:ni),stat=jerr)                                         3212
      if(jerr.ne.0) return                                                 3213
      allocate(w(1:no),stat=jerr)                                          3214
      if(jerr.ne.0) return                                                 3215
      allocate(qy(1:no),stat=jerr)                                         3216
      if(jerr.ne.0) return                                                 3217
      bta=parm                                                             3217
      omb=1.0-bta                                                          3217
      fmax=log(huge(bta)*0.1)                                              3218
      qy=q*y                                                               3218
      yb=sum(qy)                                                           3219
      if(nonzero(no,g) .ne. 0)goto 18811                                   3219
      t=0.0                                                                3220
      if(intr .eq. 0)goto 18831                                            3220
      w=q*yb                                                               3220
      az=log(yb)                                                           3220
      uu=az                                                                3221
      xm=yb*xb                                                             3221
      dv0=yb*(az-1.0)                                                      3222
      goto 18841                                                           3223
18831 continue                                                             3223
      w=q                                                                  3223
      xm=0.0                                                               3223
      uu=0.0                                                               3223
      az=uu                                                                3223
      dv0=-1.0                                                             3223
18841 continue                                                             3224
18821 continue                                                             3224
      goto 18851                                                           3225
18811 continue                                                             3225
      w=q*exp(sign(min(abs(g),fmax),g))                                    3225
      ww=sum(w)                                                            3225
      t=g                                                                  3226
      if(intr .eq. 0)goto 18871                                            3226
      eaz=yb/ww                                                            3227
      w=eaz*w                                                              3227
      az=log(eaz)                                                          3227
      uu=az                                                                3227
      dv0=dot_product(qy,g)-yb*(1.0-az)                                    3228
      goto 18881                                                           3229
18871 continue                                                             3229
      uu=0.0                                                               3229
      az=uu                                                                3229
      dv0=dot_product(qy,g)-ww                                             3229
18881 continue                                                             3230
18861 continue                                                             3230
18890 do 18891 j=1,ni                                                      3230
      if(ju(j).eq.0)goto 18891                                             3230
      jb=ix(j)                                                             3230
      je=ix(j+1)-1                                                         3231
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             3232
18891 continue                                                             3233
18892 continue                                                             3233
18851 continue                                                             3234
18801 continue                                                             3234
      tt=yb*uu                                                             3234
      ww=1.0                                                               3234
      if(intr.ne.0) ww=yb                                                  3234
      wr=qy-q*(yb*(1.0-uu))                                                3234
      a=0.0                                                                3234
      as=0.0                                                               3235
      dvr=-yb                                                              3236
18900 do 18901 i=1,no                                                      3236
      if(qy(i).gt.0.0) dvr=dvr+qy(i)*log(y(i))                             3236
18901 continue                                                             3236
18902 continue                                                             3236
      dvr=dvr-dv0                                                          3236
      dev0=dvr                                                             3238
      alf=1.0                                                              3240
      if(flmin .ge. 1.0)goto 18921                                         3240
      eqs=max(eps,flmin)                                                   3240
      alf=eqs**(1.0/(nlam-1))                                              3240
18921 continue                                                             3241
      m=0                                                                  3241
      mm=0                                                                 3241
      nlp=0                                                                3241
      nin=nlp                                                              3241
      mnl=min(mnlam,nlam)                                                  3241
      shr=shri*dev0                                                        3241
      al=0.0                                                               3241
      ixx=0                                                                3242
18930 do 18931 j=1,ni                                                      3242
      if(ju(j).eq.0)goto 18931                                             3243
      jb=ix(j)                                                             3243
      je=ix(j+1)-1                                                         3244
      ga(j)=abs(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(j)-ww*xb(j)   3246 
     *)-xb(j)*tt)/xs(j)
18931 continue                                                             3247
18932 continue                                                             3247
18940 do 18941 ilm=1,nlam                                                  3247
      al0=al                                                               3248
      if(flmin .lt. 1.0)goto 18961                                         3248
      al=ulam(ilm)                                                         3248
      goto 18951                                                           3249
18961 if(ilm .le. 2)goto 18971                                             3249
      al=al*alf                                                            3249
      goto 18951                                                           3250
18971 if(ilm .ne. 1)goto 18981                                             3250
      al=big                                                               3250
      goto 18991                                                           3251
18981 continue                                                             3251
      al0=0.0                                                              3252
19000 do 19001 j=1,ni                                                      3252
      if(ju(j).eq.0)goto 19001                                             3252
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            3252
19001 continue                                                             3253
19002 continue                                                             3253
      al0=al0/max(bta,1.0d-3)                                              3253
      al=alf*al0                                                           3254
18991 continue                                                             3255
18951 continue                                                             3255
      al2=al*omb                                                           3255
      al1=al*bta                                                           3255
      tlam=bta*(2.0*al-al0)                                                3256
19010 do 19011 k=1,ni                                                      3256
      if(ixx(k).eq.1)goto 19011                                            3256
      if(ju(k).eq.0)goto 19011                                             3257
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     3258
19011 continue                                                             3259
19012 continue                                                             3259
11020 continue                                                             3260
19020 continue                                                             3260
19021 continue                                                             3260
      az0=az                                                               3261
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                3262
19030 do 19031 j=1,ni                                                      3262
      if(ixx(j).eq.0)goto 19031                                            3262
      jb=ix(j)                                                             3262
      je=ix(j+1)-1                                                         3263
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             3264
      v(j)=(dot_product(w(jx(jb:je)),x(jb:je)**2)  -2.0*xb(j)*xm(j)+ww*x   3266 
     *b(j)**2)/xs(j)**2
19031 continue                                                             3267
19032 continue                                                             3267
19040 continue                                                             3267
19041 continue                                                             3267
      nlp=nlp+1                                                            3268
      dlx=0.0                                                              3269
19050 do 19051 k=1,ni                                                      3269
      if(ixx(k).eq.0)goto 19051                                            3269
      jb=ix(k)                                                             3269
      je=ix(k+1)-1                                                         3269
      ak=a(k)                                                              3270
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)   3272 
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                                  3273
      if(au .gt. 0.0)goto 19071                                            3273
      a(k)=0.0                                                             3273
      goto 19081                                                           3274
19071 continue                                                             3275
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))           3276
19081 continue                                                             3277
19061 continue                                                             3277
      if(a(k).eq.ak)goto 19051                                             3278
      if(mm(k) .ne. 0)goto 19101                                           3278
      nin=nin+1                                                            3278
      if(nin.gt.nx)goto 19052                                              3279
      mm(k)=nin                                                            3279
      m(nin)=k                                                             3280
19101 continue                                                             3281
      d=a(k)-ak                                                            3281
      dlx=max(dlx,v(k)*d**2)                                               3281
      dv=d/xs(k)                                                           3282
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)                 3283
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                                3284
      uu=uu-dv*xb(k)                                                       3284
      tt=tt-dv*xm(k)                                                       3285
19051 continue                                                             3286
19052 continue                                                             3286
      if(nin.gt.nx)goto 19042                                              3287
      if(intr .eq. 0)goto 19121                                            3287
      d=tt/ww-uu                                                           3288
      az=az+d                                                              3288
      dlx=max(dlx,ww*d**2)                                                 3288
      uu=uu+d                                                              3289
19121 continue                                                             3290
      if(dlx.lt.shr)goto 19042                                             3290
      if(nlp .le. maxit)goto 19141                                         3290
      jerr=-ilm                                                            3290
      return                                                               3290
19141 continue                                                             3291
19150 continue                                                             3291
19151 continue                                                             3291
      nlp=nlp+1                                                            3291
      dlx=0.0                                                              3292
19160 do 19161 l=1,nin                                                     3292
      k=m(l)                                                               3293
      jb=ix(k)                                                             3293
      je=ix(k+1)-1                                                         3293
      ak=a(k)                                                              3294
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)   3296 
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                                  3297
      if(au .gt. 0.0)goto 19181                                            3297
      a(k)=0.0                                                             3297
      goto 19191                                                           3298
19181 continue                                                             3299
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))           3300
19191 continue                                                             3301
19171 continue                                                             3301
      if(a(k).eq.ak)goto 19161                                             3301
      d=a(k)-ak                                                            3301
      dlx=max(dlx,v(k)*d**2)                                               3302
      dv=d/xs(k)                                                           3302
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)                 3303
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                                3304
      uu=uu-dv*xb(k)                                                       3304
      tt=tt-dv*xm(k)                                                       3305
19161 continue                                                             3306
19162 continue                                                             3306
      if(intr .eq. 0)goto 19211                                            3306
      d=tt/ww-uu                                                           3306
      az=az+d                                                              3307
      dlx=max(dlx,ww*d**2)                                                 3307
      uu=uu+d                                                              3308
19211 continue                                                             3309
      if(dlx.lt.shr)goto 19152                                             3309
      if(nlp .le. maxit)goto 19231                                         3309
      jerr=-ilm                                                            3309
      return                                                               3309
19231 continue                                                             3310
      goto 19151                                                           3311
19152 continue                                                             3311
      goto 19041                                                           3312
19042 continue                                                             3312
      if(nin.gt.nx)goto 19022                                              3313
      euu=exp(sign(min(abs(uu),fmax),uu))                                  3314
      w=euu*q*exp(sign(min(abs(t),fmax),t))                                3314
      ww=sum(w)                                                            3315
      wr=qy-w*(1.0-uu)                                                     3315
      tt=sum(wr)                                                           3316
      if(ww*(az-az0)**2 .ge. shr)goto 19251                                3316
      kx=0                                                                 3317
19260 do 19261 j=1,nin                                                     3317
      k=m(j)                                                               3318
      if(v(k)*(a(k)-as(k))**2.lt.shr)goto 19261                            3318
      kx=1                                                                 3318
      goto 19262                                                           3319
19261 continue                                                             3320
19262 continue                                                             3320
      if(kx .ne. 0)goto 19281                                              3321
19290 do 19291 j=1,ni                                                      3321
      if(ixx(j).eq.1)goto 19291                                            3321
      if(ju(j).eq.0)goto 19291                                             3322
      jb=ix(j)                                                             3322
      je=ix(j+1)-1                                                         3323
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             3324
      ga(j)=abs(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(j)-ww*xb(j)   3326 
     *)-xb(j)*tt)/xs(j)
      if(ga(j) .le. al1*vp(j))goto 19311                                   3326
      ixx(j)=1                                                             3326
      kx=1                                                                 3326
19311 continue                                                             3327
19291 continue                                                             3328
19292 continue                                                             3328
      if(kx.eq.1) go to 11020                                              3329
      goto 19022                                                           3330
19281 continue                                                             3331
19251 continue                                                             3332
      goto 19021                                                           3333
19022 continue                                                             3333
      if(nin .le. nx)goto 19331                                            3333
      jerr=-10000-ilm                                                      3333
      goto 18942                                                           3333
19331 continue                                                             3334
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                               3334
      kin(ilm)=nin                                                         3335
      a0(ilm)=az                                                           3335
      alm(ilm)=al                                                          3335
      lmu=ilm                                                              3336
      dev(ilm)=(dot_product(qy,t)+yb*uu-ww-dv0)/dvr                        3337
      if(ilm.lt.mnl)goto 18941                                             3337
      if(flmin.ge.1.0)goto 18941                                           3338
      me=0                                                                 3338
19340 do 19341 j=1,nin                                                     3338
      if(ca(j,ilm).ne.0.0) me=me+1                                         3338
19341 continue                                                             3338
19342 continue                                                             3338
      if(me.gt.ne)goto 18942                                               3339
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 18942              3340
      if(dev(ilm).gt.devmax)goto 18942                                     3341
18941 continue                                                             3342
18942 continue                                                             3342
      g=t+uu                                                               3343
12320 continue                                                             3343
      deallocate(t,w,wr,v,a,qy,xm,as,mm,ga,ixx)                            3344
      return                                                               3345
      end                                                                  3346
      subroutine spdeviance(no,ni,x,ix,jx,y,g,q,nlam,a0,a,flog,jerr)       3347
      implicit double precision(a-h,o-z)                                   3348
      double precision x(*),y(no),g(no),q(no),a(ni,nlam),a0(nlam),flog(n   3349 
     *lam)
      integer ix(*),jx(*)                                                  3350
      double precision, dimension (:), allocatable :: w,f                       
      if(minval(y) .ge. 0.0)goto 19361                                     3353
      jerr=8888                                                            3353
      return                                                               3353
19361 continue                                                             3354
      allocate(w(1:no),stat=jerr)                                          3355
      if(jerr.ne.0) return                                                 3356
      allocate(f(1:no),stat=jerr)                                          3357
      if(jerr.ne.0) return                                                 3358
      w=max(0d0,q)                                                         3358
      sw=sum(w)                                                            3358
      if(sw .gt. 0.0)goto 19381                                            3358
      jerr=9999                                                            3358
      go to 12320                                                          3358
19381 continue                                                             3359
      yb=dot_product(w,y)/sw                                               3359
      fmax=log(huge(y(1))*0.1)                                             3360
19390 do 19391 lam=1,nlam                                                  3360
      f=a0(lam)                                                            3361
19400 do 19401 j=1,ni                                                      3361
      if(a(j,lam).eq.0.0)goto 19401                                        3361
      jb=ix(j)                                                             3361
      je=ix(j+1)-1                                                         3362
      f(jx(jb:je))=f(jx(jb:je))+a(j,lam)*x(jb:je)                          3363
19401 continue                                                             3364
19402 continue                                                             3364
      f=f+g                                                                3365
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                   3366
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                3367
19391 continue                                                             3368
19392 continue                                                             3368
12320 continue                                                             3368
      deallocate(w,f)                                                      3369
      return                                                               3370
      end                                                                  3371
      subroutine cspdeviance(no,x,ix,jx,y,g,q,nx,nlam,a0,ca,ia,nin,flog,   3372 
     *jerr)
      implicit double precision(a-h,o-z)                                   3373
      double precision x(*),y(no),g(no),q(no),ca(nx,nlam),a0(nlam),flog(   3374 
     *nlam)
      integer ix(*),jx(*),nin(nlam),ia(nx)                                 3375
      double precision, dimension (:), allocatable :: w,f                       
      if(minval(y) .ge. 0.0)goto 19421                                     3378
      jerr=8888                                                            3378
      return                                                               3378
19421 continue                                                             3379
      allocate(w(1:no),stat=jerr)                                          3380
      if(jerr.ne.0) return                                                 3381
      allocate(f(1:no),stat=jerr)                                          3382
      if(jerr.ne.0) return                                                 3383
      w=max(0d0,q)                                                         3383
      sw=sum(w)                                                            3383
      if(sw .gt. 0.0)goto 19441                                            3383
      jerr=9999                                                            3383
      go to 12320                                                          3383
19441 continue                                                             3384
      yb=dot_product(w,y)/sw                                               3384
      fmax=log(huge(y(1))*0.1)                                             3385
19450 do 19451 lam=1,nlam                                                  3385
      f=a0(lam)                                                            3386
19460 do 19461 k=1,nin(lam)                                                3386
      j=ia(k)                                                              3386
      jb=ix(j)                                                             3386
      je=ix(j+1)-1                                                         3387
      f(jx(jb:je))=f(jx(jb:je))+ca(k,lam)*x(jb:je)                         3388
19461 continue                                                             3389
19462 continue                                                             3389
      f=f+g                                                                3390
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                   3391
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                3392
19451 continue                                                             3393
19452 continue                                                             3393
12320 continue                                                             3393
      deallocate(w,f)                                                      3394
      return                                                               3395
      end                                                                  3396
      subroutine multelnet  (parm,no,ni,nr,x,y,w,jd,vp,cl,ne,nx,nlam,flm   3399 
     *in,ulam,thr,isd,jsd,intr,maxit,  beta0,lmu,a0,ca,ia,nin,rsq,alm,nl
     *p,jerr)
      implicit double precision(a-h,o-z)                                   3400
      double precision x(no,ni),y(no,nr),w(no),vp(ni),ca(nx,nr,nlam),bet   3401 
     *a0(ni,nr)
      double precision ulam(nlam),a0(nr,nlam),rsq(nlam),alm(nlam),cl(2,n   3402 
     *i)
      integer jd(*),ia(nx),nin(nlam)                                       3403
      double precision, dimension (:), allocatable :: vq;                       
      if(maxval(vp) .gt. 0.0)goto 19481                                    3406
      jerr=10000                                                           3406
      return                                                               3406
19481 continue                                                             3407
      allocate(vq(1:ni),stat=jerr)                                         3407
      if(jerr.ne.0) return                                                 3408
      vq=max(0d0,vp)                                                       3408
      vq=vq*ni/sum(vq)                                                     3409
      call multelnetn(parm,no,ni,nr,x,y,w,jd,vq,cl,ne,nx,nlam,flmin,ulam   3411 
     *,thr,isd,  jsd,intr,maxit,beta0,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      deallocate(vq)                                                       3412
      return                                                               3413
      end                                                                  3414
      subroutine multelnetn (parm,no,ni,nr,x,y,w,jd,vp,cl,ne,nx,nlam,flm   3416 
     *in,ulam,thr,  isd,jsd,intr,maxit,beta0,lmu,a0,ca,ia,nin,rsq,alm,nl
     *p,jerr)
      implicit double precision(a-h,o-z)                                   3417
      double precision vp(ni),x(no,ni),y(no,nr),w(no),ulam(nlam),cl(2,ni   3418 
     *)
      double precision beta0(ni,nr),ca(nx,nr,nlam),a0(nr,nlam),rsq(nlam)   3419 
     *,alm(nlam)
      integer jd(*),ia(nx),nin(nlam)                                       3420
      double precision, dimension (:), allocatable :: xm,xs,xv,ym,ys            
      integer, dimension (:), allocatable :: ju                                 
      double precision, dimension (:,:,:), allocatable :: clt                   
      allocate(clt(1:2,1:nr,1:ni),stat=jerr);                                   
      if(jerr.ne.0) return                                                      
      allocate(xm(1:ni),stat=jerr)                                         3428
      if(jerr.ne.0) return                                                 3429
      allocate(xs(1:ni),stat=jerr)                                         3430
      if(jerr.ne.0) return                                                 3431
      allocate(ym(1:nr),stat=jerr)                                         3432
      if(jerr.ne.0) return                                                 3433
      allocate(ys(1:nr),stat=jerr)                                         3434
      if(jerr.ne.0) return                                                 3435
      allocate(ju(1:ni),stat=jerr)                                         3436
      if(jerr.ne.0) return                                                 3437
      allocate(xv(1:ni),stat=jerr)                                         3438
      if(jerr.ne.0) return                                                 3439
      call chkvars(no,ni,x,ju)                                             3440
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 3441
      if(maxval(ju) .gt. 0)goto 19501                                      3441
      jerr=7777                                                            3441
      return                                                               3441
19501 continue                                                             3442
      call multstandard1(no,ni,nr,x,y,w,isd,jsd,intr,ju,xm,xs,ym,ys,xv,y   3443 
     *s0,jerr)
      if(jerr.ne.0) return                                                 3444
19510 do 19511 j=1,ni                                                      3444
19520 do 19521 k=1,nr                                                      3444
19530 do 19531 i=1,2                                                       3444
      clt(i,k,j)=cl(i,j)                                                   3444
19531 continue                                                             3444
19532 continue                                                             3444
19521 continue                                                             3444
19522 continue                                                             3444
19511 continue                                                             3445
19512 continue                                                             3445
      if(isd .le. 0)goto 19551                                             3445
19560 do 19561 j=1,ni                                                      3445
19570 do 19571 k=1,nr                                                      3445
19580 do 19581 i=1,2                                                       3445
      clt(i,k,j)=clt(i,k,j)*xs(j)                                          3445
19581 continue                                                             3446
19582 continue                                                             3446
      beta0(j,k)=beta0(j,k)*xs(j)                                          3446
19571 continue                                                             3446
19572 continue                                                             3446
19561 continue                                                             3446
19562 continue                                                             3446
19551 continue                                                             3447
      if(jsd .le. 0)goto 19601                                             3447
19610 do 19611 j=1,ni                                                      3447
19620 do 19621 k=1,nr                                                      3447
19630 do 19631 i=1,2                                                       3447
      clt(i,k,j)=clt(i,k,j)/ys(k)                                          3447
19631 continue                                                             3448
19632 continue                                                             3448
      beta0(j,k)=beta0(j,k)/ys(k)                                          3448
19621 continue                                                             3448
19622 continue                                                             3448
19611 continue                                                             3448
19612 continue                                                             3448
19601 continue                                                             3449
      call multelnet2(parm,ni,nr,ju,vp,clt,y,no,ne,nx,x,nlam,flmin,ulam,   3451 
     *thr,maxit,  beta0,xv,ys0,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 3452
19640 do 19641 k=1,lmu                                                     3452
      nk=nin(k)                                                            3453
19650 do 19651 j=1,nr                                                      3454
19660 do 19661 l=1,nk                                                      3454
      ca(l,j,k)=ys(j)*ca(l,j,k)/xs(ia(l))                                  3454
19661 continue                                                             3455
19662 continue                                                             3455
      if(intr .ne. 0)goto 19681                                            3455
      a0(j,k)=0.0                                                          3455
      goto 19691                                                           3456
19681 continue                                                             3456
      a0(j,k)=ym(j)-dot_product(ca(1:nk,j,k),xm(ia(1:nk)))                 3456
19691 continue                                                             3457
19671 continue                                                             3457
19651 continue                                                             3458
19652 continue                                                             3458
19641 continue                                                             3459
19642 continue                                                             3459
      deallocate(xm,xs,ym,ys,ju,xv,clt)                                    3460
      return                                                               3461
      end                                                                  3462
      subroutine multstandard1  (no,ni,nr,x,y,w,isd,jsd,intr,ju,xm,xs,ym   3464 
     *,ys,xv,ys0,jerr)
      implicit double precision(a-h,o-z)                                   3465
      double precision x(no,ni),y(no,nr),w(no),xm(ni),xs(ni),xv(ni),ym(n   3466 
     *r),ys(nr)
      integer ju(ni)                                                       3467
      double precision, dimension (:), allocatable :: v                         
      allocate(v(1:no),stat=jerr)                                          3470
      if(jerr.ne.0) return                                                 3471
      w=w/sum(w)                                                           3471
      v=sqrt(w)                                                            3472
      if(intr .ne. 0)goto 19711                                            3473
19720 do 19721 j=1,ni                                                      3473
      if(ju(j).eq.0)goto 19721                                             3473
      xm(j)=0.0                                                            3473
      x(:,j)=v*x(:,j)                                                      3474
      z=dot_product(x(:,j),x(:,j))                                         3475
      if(isd .le. 0)goto 19741                                             3475
      xbq=dot_product(v,x(:,j))**2                                         3475
      vc=z-xbq                                                             3476
      xs(j)=sqrt(vc)                                                       3476
      x(:,j)=x(:,j)/xs(j)                                                  3476
      xv(j)=1.0+xbq/vc                                                     3477
      goto 19751                                                           3478
19741 continue                                                             3478
      xs(j)=1.0                                                            3478
      xv(j)=z                                                              3478
19751 continue                                                             3479
19731 continue                                                             3479
19721 continue                                                             3480
19722 continue                                                             3480
      ys0=0.0                                                              3481
19760 do 19761 j=1,nr                                                      3481
      ym(j)=0.0                                                            3481
      y(:,j)=v*y(:,j)                                                      3482
      z=dot_product(y(:,j),y(:,j))                                         3483
      if(jsd .le. 0)goto 19781                                             3483
      u=z-dot_product(v,y(:,j))**2                                         3483
      ys0=ys0+z/u                                                          3484
      ys(j)=sqrt(u)                                                        3484
      y(:,j)=y(:,j)/ys(j)                                                  3485
      goto 19791                                                           3486
19781 continue                                                             3486
      ys(j)=1.0                                                            3486
      ys0=ys0+z                                                            3486
19791 continue                                                             3487
19771 continue                                                             3487
19761 continue                                                             3488
19762 continue                                                             3488
      go to 10720                                                          3489
19711 continue                                                             3490
19800 do 19801 j=1,ni                                                      3490
      if(ju(j).eq.0)goto 19801                                             3491
      xm(j)=dot_product(w,x(:,j))                                          3491
      x(:,j)=v*(x(:,j)-xm(j))                                              3492
      xv(j)=dot_product(x(:,j),x(:,j))                                     3492
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       3493
19801 continue                                                             3494
19802 continue                                                             3494
      if(isd .ne. 0)goto 19821                                             3494
      xs=1.0                                                               3494
      goto 19831                                                           3495
19821 continue                                                             3495
19840 do 19841 j=1,ni                                                      3495
      if(ju(j).eq.0)goto 19841                                             3495
      x(:,j)=x(:,j)/xs(j)                                                  3495
19841 continue                                                             3496
19842 continue                                                             3496
      xv=1.0                                                               3497
19831 continue                                                             3498
19811 continue                                                             3498
      ys0=0.0                                                              3499
19850 do 19851 j=1,nr                                                      3500
      ym(j)=dot_product(w,y(:,j))                                          3500
      y(:,j)=v*(y(:,j)-ym(j))                                              3501
      z=dot_product(y(:,j),y(:,j))                                         3502
      if(jsd .le. 0)goto 19871                                             3502
      ys(j)=sqrt(z)                                                        3502
      y(:,j)=y(:,j)/ys(j)                                                  3502
      goto 19881                                                           3503
19871 continue                                                             3503
      ys0=ys0+z                                                            3503
19881 continue                                                             3504
19861 continue                                                             3504
19851 continue                                                             3505
19852 continue                                                             3505
      if(jsd .ne. 0)goto 19901                                             3505
      ys=1.0                                                               3505
      goto 19911                                                           3505
19901 continue                                                             3505
      ys0=nr                                                               3505
19911 continue                                                             3506
19891 continue                                                             3506
10720 continue                                                             3506
      deallocate(v)                                                        3507
      return                                                               3508
      end                                                                  3509
      subroutine multelnet2(beta,ni,nr,ju,vp,cl,y,no,ne,nx,x,nlam,flmin,   3511 
     *ulam,thri,  maxit,beta0,xv,ys0,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      implicit double precision(a-h,o-z)                                   3512
      double precision vp(ni),y(no,nr),x(no,ni),ulam(nlam),ao(nx,nr,nlam   3513 
     *)
      double precision beta0(ni,nr)                                        3514
      double precision rsqo(nlam),almo(nlam),xv(ni),cl(2,nr,ni)            3515
      integer ju(ni),ia(nx),kin(nlam)                                      3516
      double precision, dimension (:), allocatable :: g,gk,del,gj               
      integer, dimension (:), allocatable :: mm,ix,isc                          
      double precision, dimension (:,:), allocatable :: a                       
      allocate(a(1:nr,1:ni),stat=jerr)                                          
      if(jerr.ne.0) return                                                      
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)               3524
      allocate(gj(1:nr),stat=jerr)                                         3525
      if(jerr.ne.0) return                                                 3526
      allocate(gk(1:nr),stat=jerr)                                         3527
      if(jerr.ne.0) return                                                 3528
      allocate(del(1:nr),stat=jerr)                                        3529
      if(jerr.ne.0) return                                                 3530
      allocate(mm(1:ni),stat=jerr)                                         3531
      if(jerr.ne.0) return                                                 3532
      allocate(g(1:ni),stat=jerr)                                          3533
      if(jerr.ne.0) return                                                 3534
      allocate(ix(1:ni),stat=jerr)                                         3535
      if(jerr.ne.0) return                                                 3536
      allocate(isc(1:nr),stat=jerr)                                        3537
      if(jerr.ne.0) return                                                 3538
      bta=beta                                                             3538
      omb=1.0-bta                                                          3538
      ix=0                                                                 3538
      thr=thri*ys0/nr                                                      3540
      alf=1.0                                                              3542
      if(flmin .ge. 1.0)goto 19931                                         3542
      eqs=max(eps,flmin)                                                   3542
      alf=eqs**(1.0/(nlam-1))                                              3542
19931 continue                                                             3543
      rsq=ys0                                                              3543
      a=0.0                                                                3543
      mm=0                                                                 3543
      nlp=0                                                                3543
      nin=nlp                                                              3543
      iz=0                                                                 3543
      mnl=min(mnlam,nlam)                                                  3543
      alm=0.0                                                              3544
      if(flmin .lt. 1.0)goto 19951                                         3544
19960 do 19961 j=1,ni                                                      3544
      if(ju(j).eq.0)goto 19961                                             3545
19970 do 19971 k=1,nr                                                      3545
      y(:,k)=y(:,k)-beta0(j,k)*x(:,j)                                      3545
19971 continue                                                             3545
19972 continue                                                             3545
19961 continue                                                             3545
19962 continue                                                             3545
19951 continue                                                             3546
19980 do 19981 j=1,ni                                                      3546
      if(ju(j).eq.0)goto 19981                                             3546
      g(j)=0.0                                                             3547
19990 do 19991 k=1,nr                                                      3547
      g(j)=g(j)+dot_product(y(:,k),x(:,j))**2                              3547
19991 continue                                                             3548
19992 continue                                                             3548
      g(j)=sqrt(g(j))                                                      3549
19981 continue                                                             3550
19982 continue                                                             3550
20000 do 20001 m=1,nlam                                                    3550
      alm0=alm                                                             3551
      if(flmin .lt. 1.0)goto 20021                                         3551
      alm=ulam(m)                                                          3551
      goto 20011                                                           3552
20021 if(m .le. 2)goto 20031                                               3552
      alm=alm*alf                                                          3552
      goto 20011                                                           3553
20031 if(m .ne. 1)goto 20041                                               3553
      alm=big                                                              3553
      goto 20051                                                           3554
20041 continue                                                             3554
      alm0=0.0                                                             3555
20060 do 20061 j=1,ni                                                      3555
      if(ju(j).eq.0)goto 20061                                             3556
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           3557
20061 continue                                                             3558
20062 continue                                                             3558
      alm0=alm0/max(bta,1.0d-3)                                            3558
      alm=alf*alm0                                                         3559
20051 continue                                                             3560
20011 continue                                                             3560
      dem=alm*omb                                                          3560
      ab=alm*bta                                                           3560
      rsq0=rsq                                                             3560
      jz=1                                                                 3561
      if(flmin .lt. 1.0 .or. m .ne. 1)goto 20081                           3562
20090 do 20091 j=1,ni                                                      3562
      if(ju(j).eq.0)goto 20091                                             3562
      inz=0                                                                3563
20100 do 20101 k=1,nr                                                      3563
      a(k,j)=beta0(j,k)                                                    3563
      if(abs(beta0(j,k)).gt.1e-24) inz=1                                   3563
20101 continue                                                             3564
20102 continue                                                             3564
      if(inz .le. 0)goto 20121                                             3564
      ix(j)=1                                                              3564
      nin=nin+1                                                            3564
      mm(j)=nin                                                            3564
      ia(nin)=j                                                            3564
20121 continue                                                             3565
20091 continue                                                             3566
20092 continue                                                             3566
      goto 20131                                                           3567
20081 continue                                                             3568
      tlam=bta*(2.0*alm-alm0)                                              3569
20140 do 20141 k=1,ni                                                      3569
      if(ix(k).eq.1)goto 20141                                             3569
      if(ju(k).eq.0)goto 20141                                             3570
      if(g(k).gt.tlam*vp(k)) ix(k)=1                                       3571
20141 continue                                                             3572
20142 continue                                                             3572
20131 continue                                                             3573
20071 continue                                                             3573
20150 continue                                                             3573
20151 continue                                                             3573
      if(iz*jz.ne.0) go to 10360                                           3574
11020 continue                                                             3574
      nlp=nlp+1                                                            3574
      dlx=0.0                                                              3575
20160 do 20161 k=1,ni                                                      3575
      if(ix(k).eq.0)goto 20161                                             3575
      gkn=0.0                                                              3576
20170 do 20171 j=1,nr                                                      3576
      gj(j)=dot_product(y(:,j),x(:,k))                                     3577
      gk(j)=gj(j)+a(j,k)*xv(k)                                             3577
      gkn=gkn+gk(j)**2                                                     3579
20171 continue                                                             3579
20172 continue                                                             3579
      gkn=sqrt(gkn)                                                        3579
      u=1.0-ab*vp(k)/gkn                                                   3579
      del=a(:,k)                                                           3580
      if(u .gt. 0.0)goto 20191                                             3580
      a(:,k)=0.0                                                           3580
      goto 20201                                                           3581
20191 continue                                                             3581
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      3582
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)   3584 
     *,isc,jerr)
      if(jerr.ne.0) return                                                 3585
20201 continue                                                             3586
20181 continue                                                             3586
      del=a(:,k)-del                                                       3586
      if(maxval(abs(del)).le.0.0)goto 20161                                3587
20210 do 20211 j=1,nr                                                      3587
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              3588
      y(:,j)=y(:,j)-del(j)*x(:,k)                                          3588
      dlx=max(dlx,xv(k)*del(j)**2)                                         3589
20211 continue                                                             3590
20212 continue                                                             3590
      if(mm(k) .ne. 0)goto 20231                                           3590
      nin=nin+1                                                            3590
      if(nin.gt.nx)goto 20162                                              3591
      mm(k)=nin                                                            3591
      ia(nin)=k                                                            3592
20231 continue                                                             3593
20161 continue                                                             3594
20162 continue                                                             3594
      if(nin.gt.nx)goto 20152                                              3595
      if(dlx .ge. thr)goto 20251                                           3595
      ixx=0                                                                3596
20260 do 20261 k=1,ni                                                      3596
      if(ix(k).eq.1)goto 20261                                             3596
      if(ju(k).eq.0)goto 20261                                             3596
      g(k)=0.0                                                             3597
20270 do 20271 j=1,nr                                                      3597
      g(k)=g(k)+dot_product(y(:,j),x(:,k))**2                              3597
20271 continue                                                             3598
20272 continue                                                             3598
      g(k)=sqrt(g(k))                                                      3599
      if(g(k) .le. ab*vp(k))goto 20291                                     3599
      ix(k)=1                                                              3599
      ixx=1                                                                3599
20291 continue                                                             3600
20261 continue                                                             3601
20262 continue                                                             3601
      if(ixx.eq.1) go to 11020                                             3602
      goto 20152                                                           3603
20251 continue                                                             3604
      if(nlp .le. maxit)goto 20311                                         3604
      jerr=-m                                                              3604
      return                                                               3604
20311 continue                                                             3605
10360 continue                                                             3605
      iz=1                                                                 3606
20320 continue                                                             3606
20321 continue                                                             3606
      nlp=nlp+1                                                            3606
      dlx=0.0                                                              3607
20330 do 20331 l=1,nin                                                     3607
      k=ia(l)                                                              3607
      gkn=0.0                                                              3608
20340 do 20341 j=1,nr                                                      3608
      gj(j)=dot_product(y(:,j),x(:,k))                                     3609
      gk(j)=gj(j)+a(j,k)*xv(k)                                             3609
      gkn=gkn+gk(j)**2                                                     3611
20341 continue                                                             3611
20342 continue                                                             3611
      gkn=sqrt(gkn)                                                        3611
      u=1.0-ab*vp(k)/gkn                                                   3611
      del=a(:,k)                                                           3612
      if(u .gt. 0.0)goto 20361                                             3612
      a(:,k)=0.0                                                           3612
      goto 20371                                                           3613
20361 continue                                                             3613
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      3614
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)   3616 
     *,isc,jerr)
      if(jerr.ne.0) return                                                 3617
20371 continue                                                             3618
20351 continue                                                             3618
      del=a(:,k)-del                                                       3618
      if(maxval(abs(del)).le.0.0)goto 20331                                3619
20380 do 20381 j=1,nr                                                      3619
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              3620
      y(:,j)=y(:,j)-del(j)*x(:,k)                                          3620
      dlx=max(dlx,xv(k)*del(j)**2)                                         3621
20381 continue                                                             3622
20382 continue                                                             3622
20331 continue                                                             3623
20332 continue                                                             3623
      if(dlx.lt.thr)goto 20322                                             3623
      if(nlp .le. maxit)goto 20401                                         3623
      jerr=-m                                                              3623
      return                                                               3623
20401 continue                                                             3624
      goto 20321                                                           3625
20322 continue                                                             3625
      jz=0                                                                 3626
      goto 20151                                                           3627
20152 continue                                                             3627
      if(nin .le. nx)goto 20421                                            3627
      jerr=-10000-m                                                        3627
      goto 20002                                                           3627
20421 continue                                                             3628
      if(nin .le. 0)goto 20441                                             3628
20450 do 20451 j=1,nr                                                      3628
      ao(1:nin,j,m)=a(j,ia(1:nin))                                         3628
20451 continue                                                             3628
20452 continue                                                             3628
20441 continue                                                             3629
      kin(m)=nin                                                           3630
      rsqo(m)=1.0-rsq/ys0                                                  3630
      almo(m)=alm                                                          3630
      lmu=m                                                                3631
      if(m.lt.mnl)goto 20001                                               3631
      if(flmin.ge.1.0)goto 20001                                           3632
      me=0                                                                 3632
20460 do 20461 j=1,nin                                                     3632
      if(ao(j,1,m).ne.0.0) me=me+1                                         3632
20461 continue                                                             3632
20462 continue                                                             3632
      if(me.gt.ne)goto 20002                                               3633
      if(rsq0-rsq.lt.sml*rsq)goto 20002                                    3633
      if(rsqo(m).gt.rsqmax)goto 20002                                      3634
20001 continue                                                             3635
20002 continue                                                             3635
      deallocate(a,mm,g,ix,del,gj,gk)                                      3636
      return                                                               3637
      end                                                                  3638
      subroutine chkbnds(nr,gk,gkn,xv,cl,al1,al2,a,isc,jerr)               3639
      implicit double precision(a-h,o-z)                                   3640
      double precision gk(nr),cl(2,nr),a(nr)                               3640
      integer isc(nr)                                                      3641
      kerr=0                                                               3641
      al1p=1.0+al1/xv                                                      3641
      al2p=al2/xv                                                          3641
      isc=0                                                                3642
      gsq=gkn**2                                                           3642
      asq=dot_product(a,a)                                                 3642
      usq=0.0                                                              3644
      u=0.0                                                                3644
      kn=-1                                                                3646
20470 continue                                                             3646
20471 continue                                                             3646
      vmx=0.0                                                              3647
20480 do 20481 k=1,nr                                                      3647
      v=max(a(k)-cl(2,k),cl(1,k)-a(k))                                     3648
      if(v .le. vmx)goto 20501                                             3648
      vmx=v                                                                3648
      kn=k                                                                 3648
20501 continue                                                             3649
20481 continue                                                             3650
20482 continue                                                             3650
      if(vmx.le.0.0)goto 20472                                             3650
      if(isc(kn).ne.0)goto 20472                                           3651
      gsq=gsq-gk(kn)**2                                                    3651
      g=sqrt(gsq)/xv                                                       3652
      if(a(kn).lt.cl(1,kn)) u=cl(1,kn)                                     3652
      if(a(kn).gt.cl(2,kn)) u=cl(2,kn)                                     3653
      usq=usq+u**2                                                         3654
      if(usq .ne. 0.0)goto 20521                                           3654
      b=max(0d0,(g-al2p)/al1p)                                             3654
      goto 20531                                                           3655
20521 continue                                                             3655
      b0=sqrt(asq-a(kn)**2)                                                3656
      b=bnorm(b0,al1p,al2p,g,usq,kerr)                                     3656
      if(kerr.ne.0)goto 20472                                              3657
20531 continue                                                             3658
20511 continue                                                             3658
      asq=usq+b**2                                                         3658
      if(asq .gt. 0.0)goto 20551                                           3658
      a=0.0                                                                3658
      goto 20472                                                           3658
20551 continue                                                             3659
      a(kn)=u                                                              3659
      isc(kn)=1                                                            3659
      f=1.0/(xv*(al1p+al2p/sqrt(asq)))                                     3660
20560 do 20561 j=1,nr                                                      3660
      if(isc(j).eq.0) a(j)=f*gk(j)                                         3660
20561 continue                                                             3661
20562 continue                                                             3661
      goto 20471                                                           3662
20472 continue                                                             3662
      if(kerr.ne.0) jerr=kerr                                              3663
      return                                                               3664
      end                                                                  3665
      subroutine chkbnds1(nr,gk,gkn,xv,cl1,cl2,al1,al2,a,isc,jerr)         3666
      implicit double precision(a-h,o-z)                                   3667
      double precision gk(nr),a(nr)                                        3667
      integer isc(nr)                                                      3668
      kerr=0                                                               3668
      al1p=1.0+al1/xv                                                      3668
      al2p=al2/xv                                                          3668
      isc=0                                                                3669
      gsq=gkn**2                                                           3669
      asq=dot_product(a,a)                                                 3669
      usq=0.0                                                              3671
      u=0.0                                                                3671
      kn=-1                                                                3673
20570 continue                                                             3673
20571 continue                                                             3673
      vmx=0.0                                                              3674
20580 do 20581 k=1,nr                                                      3674
      v=max(a(k)-cl2,cl1-a(k))                                             3675
      if(v .le. vmx)goto 20601                                             3675
      vmx=v                                                                3675
      kn=k                                                                 3675
20601 continue                                                             3676
20581 continue                                                             3677
20582 continue                                                             3677
      if(vmx.le.0.0)goto 20572                                             3677
      if(isc(kn).ne.0)goto 20572                                           3678
      gsq=gsq-gk(kn)**2                                                    3678
      g=sqrt(gsq)/xv                                                       3679
      if(a(kn).lt.cl1) u=cl1                                               3679
      if(a(kn).gt.cl2) u=cl2                                               3680
      usq=usq+u**2                                                         3681
      if(usq .ne. 0.0)goto 20621                                           3681
      b=max(0d0,(g-al2p)/al1p)                                             3681
      goto 20631                                                           3682
20621 continue                                                             3682
      b0=sqrt(asq-a(kn)**2)                                                3683
      b=bnorm(b0,al1p,al2p,g,usq,kerr)                                     3683
      if(kerr.ne.0)goto 20572                                              3684
20631 continue                                                             3685
20611 continue                                                             3685
      asq=usq+b**2                                                         3685
      if(asq .gt. 0.0)goto 20651                                           3685
      a=0.0                                                                3685
      goto 20572                                                           3685
20651 continue                                                             3686
      a(kn)=u                                                              3686
      isc(kn)=1                                                            3686
      f=1.0/(xv*(al1p+al2p/sqrt(asq)))                                     3687
20660 do 20661 j=1,nr                                                      3687
      if(isc(j).eq.0) a(j)=f*gk(j)                                         3687
20661 continue                                                             3688
20662 continue                                                             3688
      goto 20571                                                           3689
20572 continue                                                             3689
      if(kerr.ne.0) jerr=kerr                                              3690
      return                                                               3691
      end                                                                  3692
      function bnorm(b0,al1p,al2p,g,usq,jerr)                              3693
      implicit double precision(a-h,o-z)                                   3694
      data thr,mxit /1.0d-10,100/                                          3695
      b=b0                                                                 3695
      zsq=b**2+usq                                                         3695
      if(zsq .gt. 0.0)goto 20681                                           3695
      bnorm=0.0                                                            3695
      return                                                               3695
20681 continue                                                             3696
      z=sqrt(zsq)                                                          3696
      f=b*(al1p+al2p/z)-g                                                  3696
      jerr=0                                                               3697
20690 do 20691 it=1,mxit                                                   3697
      b=b-f/(al1p+al2p*usq/(z*zsq))                                        3698
      zsq=b**2+usq                                                         3698
      if(zsq .gt. 0.0)goto 20711                                           3698
      bnorm=0.0                                                            3698
      return                                                               3698
20711 continue                                                             3699
      z=sqrt(zsq)                                                          3699
      f=b*(al1p+al2p/z)-g                                                  3700
      if(abs(f).le.thr)goto 20692                                          3700
      if(b .gt. 0.0)goto 20731                                             3700
      b=0.0                                                                3700
      goto 20692                                                           3700
20731 continue                                                             3701
20691 continue                                                             3702
20692 continue                                                             3702
      bnorm=b                                                              3702
      if(it.ge.mxit) jerr=90000                                            3703
      return                                                               3705
      entry chg_bnorm(arg,irg)                                             3705
      bnorm = 0.0                                                          3705
      thr=arg                                                              3705
      mxit=irg                                                             3705
      return                                                               3706
      entry get_bnorm(arg,irg)                                             3706
      bnorm = 0.0                                                          3706
      arg=thr                                                              3706
      irg=mxit                                                             3706
      return                                                               3708
      end                                                                  3709
      subroutine multsolns(ni,nx,nr,lmu,a,ia,nin,b)                        3710
      implicit double precision(a-h,o-z)                                   3711
      double precision a(nx,nr,lmu),b(ni,nr,lmu)                           3711
      integer ia(nx),nin(lmu)                                              3712
20740 do 20741 lam=1,lmu                                                   3712
      call multuncomp(ni,nr,nx,a(1,1,lam),ia,nin(lam),b(1,1,lam))          3712
20741 continue                                                             3713
20742 continue                                                             3713
      return                                                               3714
      end                                                                  3715
      subroutine multuncomp(ni,nr,nx,ca,ia,nin,a)                          3716
      implicit double precision(a-h,o-z)                                   3717
      double precision ca(nx,nr),a(ni,nr)                                  3717
      integer ia(nx)                                                       3718
      a=0.0                                                                3719
      if(nin .le. 0)goto 20761                                             3719
20770 do 20771 j=1,nr                                                      3719
      a(ia(1:nin),j)=ca(1:nin,j)                                           3719
20771 continue                                                             3719
20772 continue                                                             3719
20761 continue                                                             3720
      return                                                               3721
      end                                                                  3722
      subroutine multmodval(nx,nr,a0,ca,ia,nin,n,x,f)                      3723
      implicit double precision(a-h,o-z)                                   3724
      double precision a0(nr),ca(nx,nr),x(n,*),f(nr,n)                     3724
      integer ia(nx)                                                       3725
20780 do 20781 i=1,n                                                       3725
      f(:,i)=a0                                                            3725
20781 continue                                                             3725
20782 continue                                                             3725
      if(nin.le.0) return                                                  3726
20790 do 20791 i=1,n                                                       3726
20800 do 20801 j=1,nr                                                      3726
      f(j,i)=f(j,i)+dot_product(ca(1:nin,j),x(i,ia(1:nin)))                3726
20801 continue                                                             3726
20802 continue                                                             3726
20791 continue                                                             3727
20792 continue                                                             3727
      return                                                               3728
      end                                                                  3729
      subroutine multspelnet  (parm,no,ni,nr,x,ix,jx,y,w,jd,vp,cl,ne,nx,   3732 
     *nlam,flmin,ulam,thr,isd,  jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,
     *nlp,jerr)
      implicit double precision(a-h,o-z)                                   3733
      double precision x(*),y(no,nr),w(no),vp(ni),ulam(nlam),cl(2,ni)      3734
      double precision ca(nx,nr,nlam),a0(nr,nlam),rsq(nlam),alm(nlam)      3735
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           3736
      double precision, dimension (:), allocatable :: vq;                       
      if(maxval(vp) .gt. 0.0)goto 20821                                    3739
      jerr=10000                                                           3739
      return                                                               3739
20821 continue                                                             3740
      allocate(vq(1:ni),stat=jerr)                                         3740
      if(jerr.ne.0) return                                                 3741
      vq=max(0d0,vp)                                                       3741
      vq=vq*ni/sum(vq)                                                     3742
      call multspelnetn(parm,no,ni,nr,x,ix,jx,y,w,jd,vq,cl,ne,nx,nlam,fl   3744 
     *min,  ulam,thr,isd,jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jer
     *r)
      deallocate(vq)                                                       3745
      return                                                               3746
      end                                                                  3747
      subroutine multspelnetn(parm,no,ni,nr,x,ix,jx,y,w,jd,vp,cl,ne,nx,n   3749 
     *lam,flmin,  ulam,thr,isd,jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,n
     *lp,jerr)
      implicit double precision(a-h,o-z)                                   3750
      double precision x(*),vp(ni),y(no,nr),w(no),ulam(nlam),cl(2,ni)      3751
      double precision ca(nx,nr,nlam),a0(nr,nlam),rsq(nlam),alm(nlam)      3752
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           3753
      double precision, dimension (:), allocatable :: xm,xs,xv,ym,ys            
      integer, dimension (:), allocatable :: ju                                 
      double precision, dimension (:,:,:), allocatable :: clt                   
      allocate(clt(1:2,1:nr,1:ni),stat=jerr)                                    
      if(jerr.ne.0) return                                                      
      allocate(xm(1:ni),stat=jerr)                                         3761
      if(jerr.ne.0) return                                                 3762
      allocate(xs(1:ni),stat=jerr)                                         3763
      if(jerr.ne.0) return                                                 3764
      allocate(ym(1:nr),stat=jerr)                                         3765
      if(jerr.ne.0) return                                                 3766
      allocate(ys(1:nr),stat=jerr)                                         3767
      if(jerr.ne.0) return                                                 3768
      allocate(ju(1:ni),stat=jerr)                                         3769
      if(jerr.ne.0) return                                                 3770
      allocate(xv(1:ni),stat=jerr)                                         3771
      if(jerr.ne.0) return                                                 3772
      call spchkvars(no,ni,x,ix,ju)                                        3773
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 3774
      if(maxval(ju) .gt. 0)goto 20841                                      3774
      jerr=7777                                                            3774
      return                                                               3774
20841 continue                                                             3775
      call multspstandard1(no,ni,nr,x,ix,jx,y,w,ju,isd,jsd,intr,  xm,xs,   3777 
     *ym,ys,xv,ys0,jerr)
      if(jerr.ne.0) return                                                 3778
20850 do 20851 j=1,ni                                                      3778
20860 do 20861 k=1,nr                                                      3778
20870 do 20871 i=1,2                                                       3778
      clt(i,k,j)=cl(i,j)                                                   3778
20871 continue                                                             3778
20872 continue                                                             3778
20861 continue                                                             3778
20862 continue                                                             3778
20851 continue                                                             3779
20852 continue                                                             3779
      if(isd .le. 0)goto 20891                                             3779
20900 do 20901 j=1,ni                                                      3779
20910 do 20911 k=1,nr                                                      3779
20920 do 20921 i=1,2                                                       3779
      clt(i,k,j)=clt(i,k,j)*xs(j)                                          3779
20921 continue                                                             3779
20922 continue                                                             3779
20911 continue                                                             3779
20912 continue                                                             3779
20901 continue                                                             3779
20902 continue                                                             3779
20891 continue                                                             3780
      if(jsd .le. 0)goto 20941                                             3780
20950 do 20951 j=1,ni                                                      3780
20960 do 20961 k=1,nr                                                      3780
20970 do 20971 i=1,2                                                       3780
      clt(i,k,j)=clt(i,k,j)/ys(k)                                          3780
20971 continue                                                             3780
20972 continue                                                             3780
20961 continue                                                             3780
20962 continue                                                             3780
20951 continue                                                             3780
20952 continue                                                             3780
20941 continue                                                             3781
      call multspelnet2(parm,ni,nr,y,w,no,ne,nx,x,ix,jx,ju,vp,clt,nlam,f   3783 
     *lmin,  ulam,thr,maxit,xm,xs,xv,ys0,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 3784
20980 do 20981 k=1,lmu                                                     3784
      nk=nin(k)                                                            3785
20990 do 20991 j=1,nr                                                      3786
21000 do 21001 l=1,nk                                                      3786
      ca(l,j,k)=ys(j)*ca(l,j,k)/xs(ia(l))                                  3786
21001 continue                                                             3787
21002 continue                                                             3787
      if(intr .ne. 0)goto 21021                                            3787
      a0(j,k)=0.0                                                          3787
      goto 21031                                                           3788
21021 continue                                                             3788
      a0(j,k)=ym(j)-dot_product(ca(1:nk,j,k),xm(ia(1:nk)))                 3788
21031 continue                                                             3789
21011 continue                                                             3789
20991 continue                                                             3790
20992 continue                                                             3790
20981 continue                                                             3791
20982 continue                                                             3791
      deallocate(xm,xs,ym,ys,ju,xv,clt)                                    3792
      return                                                               3793
      end                                                                  3794
      subroutine multspstandard1(no,ni,nr,x,ix,jx,y,w,ju,isd,jsd,intr,     3796 
     *xm,xs,ym,ys,xv,ys0,jerr)
      implicit double precision(a-h,o-z)                                   3797
      double precision x(*),y(no,nr),w(no),xm(ni),xs(ni),xv(ni),ym(nr),y   3798 
     *s(nr)
      integer ix(*),jx(*),ju(ni)                                           3799
      w=w/sum(w)                                                           3800
      if(intr .ne. 0)goto 21051                                            3801
21060 do 21061 j=1,ni                                                      3801
      if(ju(j).eq.0)goto 21061                                             3801
      xm(j)=0.0                                                            3801
      jb=ix(j)                                                             3801
      je=ix(j+1)-1                                                         3802
      z=dot_product(w(jx(jb:je)),x(jb:je)**2)                              3803
      if(isd .le. 0)goto 21081                                             3803
      xbq=dot_product(w(jx(jb:je)),x(jb:je))**2                            3803
      vc=z-xbq                                                             3804
      xs(j)=sqrt(vc)                                                       3804
      xv(j)=1.0+xbq/vc                                                     3805
      goto 21091                                                           3806
21081 continue                                                             3806
      xs(j)=1.0                                                            3806
      xv(j)=z                                                              3806
21091 continue                                                             3807
21071 continue                                                             3807
21061 continue                                                             3808
21062 continue                                                             3808
      ys0=0.0                                                              3809
21100 do 21101 j=1,nr                                                      3809
      ym(j)=0.0                                                            3809
      z=dot_product(w,y(:,j)**2)                                           3810
      if(jsd .le. 0)goto 21121                                             3810
      u=z-dot_product(w,y(:,j))**2                                         3810
      ys0=ys0+z/u                                                          3811
      ys(j)=sqrt(u)                                                        3811
      y(:,j)=y(:,j)/ys(j)                                                  3812
      goto 21131                                                           3813
21121 continue                                                             3813
      ys(j)=1.0                                                            3813
      ys0=ys0+z                                                            3813
21131 continue                                                             3814
21111 continue                                                             3814
21101 continue                                                             3815
21102 continue                                                             3815
      return                                                               3816
21051 continue                                                             3817
21140 do 21141 j=1,ni                                                      3817
      if(ju(j).eq.0)goto 21141                                             3818
      jb=ix(j)                                                             3818
      je=ix(j+1)-1                                                         3818
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             3819
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 3820
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       3821
21141 continue                                                             3822
21142 continue                                                             3822
      if(isd .ne. 0)goto 21161                                             3822
      xs=1.0                                                               3822
      goto 21171                                                           3822
21161 continue                                                             3822
      xv=1.0                                                               3822
21171 continue                                                             3823
21151 continue                                                             3823
      ys0=0.0                                                              3824
21180 do 21181 j=1,nr                                                      3825
      ym(j)=dot_product(w,y(:,j))                                          3825
      y(:,j)=y(:,j)-ym(j)                                                  3826
      z=dot_product(w,y(:,j)**2)                                           3827
      if(jsd .le. 0)goto 21201                                             3827
      ys(j)=sqrt(z)                                                        3827
      y(:,j)=y(:,j)/ys(j)                                                  3827
      goto 21211                                                           3828
21201 continue                                                             3828
      ys0=ys0+z                                                            3828
21211 continue                                                             3829
21191 continue                                                             3829
21181 continue                                                             3830
21182 continue                                                             3830
      if(jsd .ne. 0)goto 21231                                             3830
      ys=1.0                                                               3830
      goto 21241                                                           3830
21231 continue                                                             3830
      ys0=nr                                                               3830
21241 continue                                                             3831
21221 continue                                                             3831
      return                                                               3832
      end                                                                  3833
      subroutine multspelnet2(beta,ni,nr,y,w,no,ne,nx,x,ix,jx,ju,vp,cl,n   3835 
     *lam,flmin,  ulam,thri,maxit,xm,xs,xv,ys0,lmu,ao,ia,kin,rsqo,almo,n
     *lp,jerr)
      implicit double precision(a-h,o-z)                                   3836
      double precision y(no,nr),w(no),x(*),vp(ni),ulam(nlam),cl(2,nr,ni)   3837
      double precision ao(nx,nr,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni   3838 
     *),xv(ni)
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                          3839
      double precision, dimension (:), allocatable :: g,gj,gk,del,o             
      integer, dimension (:), allocatable :: mm,iy,isc                          
      double precision, dimension (:,:), allocatable :: a                       
      allocate(a(1:nr,1:ni),stat=jerr)                                          
      if(jerr.ne.0) return                                                      
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)               3847
      allocate(mm(1:ni),stat=jerr)                                         3848
      if(jerr.ne.0) return                                                 3849
      allocate(g(1:ni),stat=jerr)                                          3850
      if(jerr.ne.0) return                                                 3851
      allocate(gj(1:nr),stat=jerr)                                         3852
      if(jerr.ne.0) return                                                 3853
      allocate(gk(1:nr),stat=jerr)                                         3854
      if(jerr.ne.0) return                                                 3855
      allocate(del(1:nr),stat=jerr)                                        3856
      if(jerr.ne.0) return                                                 3857
      allocate(o(1:nr),stat=jerr)                                          3858
      if(jerr.ne.0) return                                                 3859
      allocate(iy(1:ni),stat=jerr)                                         3860
      if(jerr.ne.0) return                                                 3861
      allocate(isc(1:nr),stat=jerr)                                        3862
      if(jerr.ne.0) return                                                 3863
      bta=beta                                                             3863
      omb=1.0-bta                                                          3863
      alm=0.0                                                              3863
      iy=0                                                                 3863
      thr=thri*ys0/nr                                                      3865
      alf=1.0                                                              3867
      if(flmin .ge. 1.0)goto 21261                                         3867
      eqs=max(eps,flmin)                                                   3867
      alf=eqs**(1.0/(nlam-1))                                              3867
21261 continue                                                             3868
      rsq=ys0                                                              3868
      a=0.0                                                                3868
      mm=0                                                                 3868
      o=0.0                                                                3868
      nlp=0                                                                3868
      nin=nlp                                                              3868
      iz=0                                                                 3868
      mnl=min(mnlam,nlam)                                                  3869
21270 do 21271 j=1,ni                                                      3869
      if(ju(j).eq.0)goto 21271                                             3869
      jb=ix(j)                                                             3869
      je=ix(j+1)-1                                                         3869
      g(j)=0.0                                                             3870
21280 do 21281 k=1,nr                                                      3871
      g(j)=g(j)+(dot_product(y(jx(jb:je),k),w(jx(jb:je))*x(jb:je))/xs(j)   3872 
     *)**2
21281 continue                                                             3873
21282 continue                                                             3873
      g(j)=sqrt(g(j))                                                      3874
21271 continue                                                             3875
21272 continue                                                             3875
21290 do 21291 m=1,nlam                                                    3875
      alm0=alm                                                             3876
      if(flmin .lt. 1.0)goto 21311                                         3876
      alm=ulam(m)                                                          3876
      goto 21301                                                           3877
21311 if(m .le. 2)goto 21321                                               3877
      alm=alm*alf                                                          3877
      goto 21301                                                           3878
21321 if(m .ne. 1)goto 21331                                               3878
      alm=big                                                              3878
      goto 21341                                                           3879
21331 continue                                                             3879
      alm0=0.0                                                             3880
21350 do 21351 j=1,ni                                                      3880
      if(ju(j).eq.0)goto 21351                                             3881
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           3882
21351 continue                                                             3883
21352 continue                                                             3883
      alm0=alm0/max(bta,1.0d-3)                                            3883
      alm=alf*alm0                                                         3884
21341 continue                                                             3885
21301 continue                                                             3885
      dem=alm*omb                                                          3885
      ab=alm*bta                                                           3885
      rsq0=rsq                                                             3885
      jz=1                                                                 3886
      tlam=bta*(2.0*alm-alm0)                                              3887
21360 do 21361 k=1,ni                                                      3887
      if(iy(k).eq.1)goto 21361                                             3887
      if(ju(k).eq.0)goto 21361                                             3888
      if(g(k).gt.tlam*vp(k)) iy(k)=1                                       3889
21361 continue                                                             3890
21362 continue                                                             3890
21370 continue                                                             3890
21371 continue                                                             3890
      if(iz*jz.ne.0) go to 10360                                           3891
11020 continue                                                             3891
      nlp=nlp+1                                                            3891
      dlx=0.0                                                              3892
21380 do 21381 k=1,ni                                                      3892
      if(iy(k).eq.0)goto 21381                                             3892
      jb=ix(k)                                                             3892
      je=ix(k+1)-1                                                         3892
      gkn=0.0                                                              3893
21390 do 21391 j=1,nr                                                      3894
      gj(j)=dot_product(y(jx(jb:je),j)+o(j),w(jx(jb:je))*x(jb:je))/xs(k)   3895
      gk(j)=gj(j)+a(j,k)*xv(k)                                             3895
      gkn=gkn+gk(j)**2                                                     3896
21391 continue                                                             3897
21392 continue                                                             3897
      gkn=sqrt(gkn)                                                        3897
      u=1.0-ab*vp(k)/gkn                                                   3897
      del=a(:,k)                                                           3898
      if(u .gt. 0.0)goto 21411                                             3898
      a(:,k)=0.0                                                           3898
      goto 21421                                                           3899
21411 continue                                                             3899
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      3900
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)   3902 
     *,isc,jerr)
      if(jerr.ne.0) return                                                 3903
21421 continue                                                             3904
21401 continue                                                             3904
      del=a(:,k)-del                                                       3904
      if(maxval(abs(del)).le.0.0)goto 21381                                3905
      if(mm(k) .ne. 0)goto 21441                                           3905
      nin=nin+1                                                            3905
      if(nin.gt.nx)goto 21382                                              3906
      mm(k)=nin                                                            3906
      ia(nin)=k                                                            3907
21441 continue                                                             3908
21450 do 21451 j=1,nr                                                      3908
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              3909
      y(jx(jb:je),j)=y(jx(jb:je),j)-del(j)*x(jb:je)/xs(k)                  3910
      o(j)=o(j)+del(j)*xm(k)/xs(k)                                         3910
      dlx=max(xv(k)*del(j)**2,dlx)                                         3911
21451 continue                                                             3912
21452 continue                                                             3912
21381 continue                                                             3913
21382 continue                                                             3913
      if(nin.gt.nx)goto 21372                                              3914
      if(dlx .ge. thr)goto 21471                                           3914
      ixx=0                                                                3915
21480 do 21481 j=1,ni                                                      3915
      if(iy(j).eq.1)goto 21481                                             3915
      if(ju(j).eq.0)goto 21481                                             3916
      jb=ix(j)                                                             3916
      je=ix(j+1)-1                                                         3916
      g(j)=0.0                                                             3917
21490 do 21491 k=1,nr                                                      3917
      g(j)=g(j)+  (dot_product(y(jx(jb:je),k)+o(k),w(jx(jb:je))*x(jb:je)   3919 
     *)/xs(j))**2
21491 continue                                                             3920
21492 continue                                                             3920
      g(j)=sqrt(g(j))                                                      3921
      if(g(j) .le. ab*vp(j))goto 21511                                     3921
      iy(j)=1                                                              3921
      ixx=1                                                                3921
21511 continue                                                             3922
21481 continue                                                             3923
21482 continue                                                             3923
      if(ixx.eq.1) go to 11020                                             3924
      goto 21372                                                           3925
21471 continue                                                             3926
      if(nlp .le. maxit)goto 21531                                         3926
      jerr=-m                                                              3926
      return                                                               3926
21531 continue                                                             3927
10360 continue                                                             3927
      iz=1                                                                 3928
21540 continue                                                             3928
21541 continue                                                             3928
      nlp=nlp+1                                                            3928
      dlx=0.0                                                              3929
21550 do 21551 l=1,nin                                                     3929
      k=ia(l)                                                              3929
      jb=ix(k)                                                             3929
      je=ix(k+1)-1                                                         3929
      gkn=0.0                                                              3930
21560 do 21561 j=1,nr                                                      3930
      gj(j)=  dot_product(y(jx(jb:je),j)+o(j),w(jx(jb:je))*x(jb:je))/xs(   3932 
     *k)
      gk(j)=gj(j)+a(j,k)*xv(k)                                             3932
      gkn=gkn+gk(j)**2                                                     3933
21561 continue                                                             3934
21562 continue                                                             3934
      gkn=sqrt(gkn)                                                        3934
      u=1.0-ab*vp(k)/gkn                                                   3934
      del=a(:,k)                                                           3935
      if(u .gt. 0.0)goto 21581                                             3935
      a(:,k)=0.0                                                           3935
      goto 21591                                                           3936
21581 continue                                                             3936
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      3937
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)   3939 
     *,isc,jerr)
      if(jerr.ne.0) return                                                 3940
21591 continue                                                             3941
21571 continue                                                             3941
      del=a(:,k)-del                                                       3941
      if(maxval(abs(del)).le.0.0)goto 21551                                3942
21600 do 21601 j=1,nr                                                      3942
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              3943
      y(jx(jb:je),j)=y(jx(jb:je),j)-del(j)*x(jb:je)/xs(k)                  3944
      o(j)=o(j)+del(j)*xm(k)/xs(k)                                         3944
      dlx=max(xv(k)*del(j)**2,dlx)                                         3945
21601 continue                                                             3946
21602 continue                                                             3946
21551 continue                                                             3947
21552 continue                                                             3947
      if(dlx.lt.thr)goto 21542                                             3947
      if(nlp .le. maxit)goto 21621                                         3947
      jerr=-m                                                              3947
      return                                                               3947
21621 continue                                                             3948
      goto 21541                                                           3949
21542 continue                                                             3949
      jz=0                                                                 3950
      goto 21371                                                           3951
21372 continue                                                             3951
      if(nin .le. nx)goto 21641                                            3951
      jerr=-10000-m                                                        3951
      goto 21292                                                           3951
21641 continue                                                             3952
      if(nin .le. 0)goto 21661                                             3952
21670 do 21671 j=1,nr                                                      3952
      ao(1:nin,j,m)=a(j,ia(1:nin))                                         3952
21671 continue                                                             3952
21672 continue                                                             3952
21661 continue                                                             3953
      kin(m)=nin                                                           3954
      rsqo(m)=1.0-rsq/ys0                                                  3954
      almo(m)=alm                                                          3954
      lmu=m                                                                3955
      if(m.lt.mnl)goto 21291                                               3955
      if(flmin.ge.1.0)goto 21291                                           3956
      me=0                                                                 3956
21680 do 21681 j=1,nin                                                     3956
      if(ao(j,1,m).ne.0.0) me=me+1                                         3956
21681 continue                                                             3956
21682 continue                                                             3956
      if(me.gt.ne)goto 21292                                               3957
      if(rsq0-rsq.lt.sml*rsq)goto 21292                                    3957
      if(rsqo(m).gt.rsqmax)goto 21292                                      3958
21291 continue                                                             3959
21292 continue                                                             3959
      deallocate(a,mm,g,iy,gj,gk,del,o)                                    3960
      return                                                               3961
      end                                                                  3962
      subroutine multlognetn(parm,no,ni,nc,x,y,g,w,ju,vp,cl,ne,nx,nlam,f   3964 
     *lmin,ulam,  shri,intr,maxit,xv,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jer
     *r)
      implicit double precision(a-h,o-z)                                   3965
      double precision x(no,ni),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam   3966 
     *),cl(2,ni)
      double precision a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),xv(   3967 
     *ni)
      integer ju(ni),m(nx),kin(nlam)                                       3968
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
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               3981
      exmn=-exmx                                                           3982
      allocate(mm(1:ni),stat=jerr)                                         3983
      if(jerr.ne.0) return                                                 3984
      allocate(is(1:max(nc,ni)),stat=jerr)                                 3985
      if(jerr.ne.0) return                                                 3986
      allocate(sxp(1:no),stat=jerr)                                        3987
      if(jerr.ne.0) return                                                 3988
      allocate(sxpl(1:no),stat=jerr)                                       3989
      if(jerr.ne.0) return                                                 3990
      allocate(ga(1:ni),stat=jerr)                                         3991
      if(jerr.ne.0) return                                                 3992
      allocate(ixx(1:ni),stat=jerr)                                        3993
      if(jerr.ne.0) return                                                 3994
      allocate(gk(1:nc),stat=jerr)                                         3995
      if(jerr.ne.0) return                                                 3996
      allocate(del(1:nc),stat=jerr)                                        3997
      if(jerr.ne.0) return                                                 3998
      allocate(isc(1:nc),stat=jerr)                                        3999
      if(jerr.ne.0) return                                                 4000
      pmax=1.0-pmin                                                        4000
      emin=pmin/pmax                                                       4000
      emax=1.0/emin                                                        4001
      bta=parm                                                             4001
      omb=1.0-bta                                                          4001
      dev1=0.0                                                             4001
      dev0=0.0                                                             4002
21690 do 21691 ic=1,nc                                                     4002
      q0=dot_product(w,y(:,ic))                                            4003
      if(q0 .gt. pmin)goto 21711                                           4003
      jerr =8000+ic                                                        4003
      return                                                               4003
21711 continue                                                             4004
      if(q0 .lt. pmax)goto 21731                                           4004
      jerr =9000+ic                                                        4004
      return                                                               4004
21731 continue                                                             4005
      if(intr .ne. 0)goto 21751                                            4005
      q0=1.0/nc                                                            4005
      b(0,ic)=0.0                                                          4005
      goto 21761                                                           4006
21751 continue                                                             4006
      b(0,ic)=log(q0)                                                      4006
      dev1=dev1-q0*b(0,ic)                                                 4006
21761 continue                                                             4007
21741 continue                                                             4007
      b(1:ni,ic)=0.0                                                       4008
21691 continue                                                             4009
21692 continue                                                             4009
      if(intr.eq.0) dev1=log(float(nc))                                    4009
      ixx=0                                                                4009
      al=0.0                                                               4010
      if(nonzero(no*nc,g) .ne. 0)goto 21781                                4011
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         4011
      sxp=0.0                                                              4012
21790 do 21791 ic=1,nc                                                     4012
      q(:,ic)=exp(b(0,ic))                                                 4012
      sxp=sxp+q(:,ic)                                                      4012
21791 continue                                                             4013
21792 continue                                                             4013
      goto 21801                                                           4014
21781 continue                                                             4014
21810 do 21811 i=1,no                                                      4014
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         4014
21811 continue                                                             4014
21812 continue                                                             4014
      sxp=0.0                                                              4015
      if(intr .ne. 0)goto 21831                                            4015
      b(0,:)=0.0                                                           4015
      goto 21841                                                           4016
21831 continue                                                             4016
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 4016
      if(jerr.ne.0) return                                                 4016
21841 continue                                                             4017
21821 continue                                                             4017
      dev1=0.0                                                             4018
21850 do 21851 ic=1,nc                                                     4018
      q(:,ic)=b(0,ic)+g(:,ic)                                              4019
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             4020
      q(:,ic)=exp(q(:,ic))                                                 4020
      sxp=sxp+q(:,ic)                                                      4021
21851 continue                                                             4022
21852 continue                                                             4022
      sxpl=w*log(sxp)                                                      4022
21860 do 21861 ic=1,nc                                                     4022
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  4022
21861 continue                                                             4023
21862 continue                                                             4023
21801 continue                                                             4024
21771 continue                                                             4024
21870 do 21871 ic=1,nc                                                     4024
21880 do 21881 i=1,no                                                      4024
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               4024
21881 continue                                                             4024
21882 continue                                                             4024
21871 continue                                                             4025
21872 continue                                                             4025
      dev0=dev0+dev1                                                       4027
      alf=1.0                                                              4029
      if(flmin .ge. 1.0)goto 21901                                         4029
      eqs=max(eps,flmin)                                                   4029
      alf=eqs**(1.0/(nlam-1))                                              4029
21901 continue                                                             4030
      m=0                                                                  4030
      mm=0                                                                 4030
      nin=0                                                                4030
      nlp=0                                                                4030
      mnl=min(mnlam,nlam)                                                  4030
      bs=0.0                                                               4030
      shr=shri*dev0                                                        4031
      ga=0.0                                                               4032
21910 do 21911 ic=1,nc                                                     4032
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      4033
21920 do 21921 j=1,ni                                                      4033
      if(ju(j).ne.0) ga(j)=ga(j)+dot_product(r(:,ic),x(:,j))**2            4033
21921 continue                                                             4034
21922 continue                                                             4034
21911 continue                                                             4035
21912 continue                                                             4035
      ga=sqrt(ga)                                                          4036
21930 do 21931 ilm=1,nlam                                                  4036
      al0=al                                                               4037
      if(flmin .lt. 1.0)goto 21951                                         4037
      al=ulam(ilm)                                                         4037
      goto 21941                                                           4038
21951 if(ilm .le. 2)goto 21961                                             4038
      al=al*alf                                                            4038
      goto 21941                                                           4039
21961 if(ilm .ne. 1)goto 21971                                             4039
      al=big                                                               4039
      goto 21981                                                           4040
21971 continue                                                             4040
      al0=0.0                                                              4041
21990 do 21991 j=1,ni                                                      4041
      if(ju(j).eq.0)goto 21991                                             4041
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            4041
21991 continue                                                             4042
21992 continue                                                             4042
      al0=al0/max(bta,1.0d-3)                                              4042
      al=alf*al0                                                           4043
21981 continue                                                             4044
21941 continue                                                             4044
      al2=al*omb                                                           4044
      al1=al*bta                                                           4044
      tlam=bta*(2.0*al-al0)                                                4045
22000 do 22001 k=1,ni                                                      4045
      if(ixx(k).eq.1)goto 22001                                            4045
      if(ju(k).eq.0)goto 22001                                             4046
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     4047
22001 continue                                                             4048
22002 continue                                                             4048
11020 continue                                                             4049
22010 continue                                                             4049
22011 continue                                                             4049
      ix=0                                                                 4049
      jx=ix                                                                4049
      kx=jx                                                                4049
      t=0.0                                                                4050
22020 do 22021 ic=1,nc                                                     4050
      t=max(t,maxval(q(:,ic)*(1.0-q(:,ic)/sxp)/sxp))                       4050
22021 continue                                                             4051
22022 continue                                                             4051
      if(t .ge. eps)goto 22041                                             4051
      kx=1                                                                 4051
      goto 22012                                                           4051
22041 continue                                                             4051
      t=2.0*t                                                              4051
      alt=al1/t                                                            4051
      al2t=al2/t                                                           4052
22050 do 22051 ic=1,nc                                                     4053
      bs(0,ic)=b(0,ic)                                                     4053
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          4054
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)/t                                    4055
      d=0.0                                                                4055
      if(intr.ne.0) d=sum(r(:,ic))                                         4056
      if(d .eq. 0.0)goto 22071                                             4057
      b(0,ic)=b(0,ic)+d                                                    4057
      r(:,ic)=r(:,ic)-d*w                                                  4057
      dlx=max(dlx,d**2)                                                    4058
22071 continue                                                             4059
22051 continue                                                             4060
22052 continue                                                             4060
22080 continue                                                             4060
22081 continue                                                             4060
      nlp=nlp+nc                                                           4060
      dlx=0.0                                                              4061
22090 do 22091 k=1,ni                                                      4061
      if(ixx(k).eq.0)goto 22091                                            4061
      gkn=0.0                                                              4062
22100 do 22101 ic=1,nc                                                     4062
      gk(ic)=dot_product(r(:,ic),x(:,k))+b(k,ic)*xv(k)                     4063
      gkn=gkn+gk(ic)**2                                                    4064
22101 continue                                                             4065
22102 continue                                                             4065
      gkn=sqrt(gkn)                                                        4065
      u=1.0-alt*vp(k)/gkn                                                  4065
      del=b(k,:)                                                           4066
      if(u .gt. 0.0)goto 22121                                             4066
      b(k,:)=0.0                                                           4066
      goto 22131                                                           4067
22121 continue                                                             4067
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                     4068
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),  cl(2,k),vp(k)*al2t,alt*vp(   4070 
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                                 4071
22131 continue                                                             4072
22111 continue                                                             4072
      del=b(k,:)-del                                                       4072
      if(maxval(abs(del)).le.0.0)goto 22091                                4073
22140 do 22141 ic=1,nc                                                     4073
      dlx=max(dlx,xv(k)*del(ic)**2)                                        4074
      r(:,ic)=r(:,ic)-del(ic)*w*x(:,k)                                     4075
22141 continue                                                             4076
22142 continue                                                             4076
      if(mm(k) .ne. 0)goto 22161                                           4076
      nin=nin+1                                                            4077
      if(nin .le. nx)goto 22181                                            4077
      jx=1                                                                 4077
      goto 22092                                                           4077
22181 continue                                                             4078
      mm(k)=nin                                                            4078
      m(nin)=k                                                             4079
22161 continue                                                             4080
22091 continue                                                             4081
22092 continue                                                             4081
      if(jx.gt.0)goto 22082                                                4081
      if(dlx.lt.shr)goto 22082                                             4082
      if(nlp .le. maxit)goto 22201                                         4082
      jerr=-ilm                                                            4082
      return                                                               4082
22201 continue                                                             4083
22210 continue                                                             4083
22211 continue                                                             4083
      nlp=nlp+nc                                                           4083
      dlx=0.0                                                              4084
22220 do 22221 l=1,nin                                                     4084
      k=m(l)                                                               4084
      gkn=0.0                                                              4085
22230 do 22231 ic=1,nc                                                     4085
      gk(ic)=dot_product(r(:,ic),x(:,k))+b(k,ic)*xv(k)                     4086
      gkn=gkn+gk(ic)**2                                                    4087
22231 continue                                                             4088
22232 continue                                                             4088
      gkn=sqrt(gkn)                                                        4088
      u=1.0-alt*vp(k)/gkn                                                  4088
      del=b(k,:)                                                           4089
      if(u .gt. 0.0)goto 22251                                             4089
      b(k,:)=0.0                                                           4089
      goto 22261                                                           4090
22251 continue                                                             4090
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                     4091
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),  cl(2,k),vp(k)*al2t,alt*vp(   4093 
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                                 4094
22261 continue                                                             4095
22241 continue                                                             4095
      del=b(k,:)-del                                                       4095
      if(maxval(abs(del)).le.0.0)goto 22221                                4096
22270 do 22271 ic=1,nc                                                     4096
      dlx=max(dlx,xv(k)*del(ic)**2)                                        4097
      r(:,ic)=r(:,ic)-del(ic)*w*x(:,k)                                     4098
22271 continue                                                             4099
22272 continue                                                             4099
22221 continue                                                             4100
22222 continue                                                             4100
      if(dlx.lt.shr)goto 22212                                             4100
      if(nlp .le. maxit)goto 22291                                         4100
      jerr=-ilm                                                            4100
      return                                                               4100
22291 continue                                                             4102
      goto 22211                                                           4103
22212 continue                                                             4103
      goto 22081                                                           4104
22082 continue                                                             4104
      if(jx.gt.0)goto 22012                                                4105
22300 do 22301 ic=1,nc                                                     4106
      if((b(0,ic)-bs(0,ic))**2.gt.shr) ix=1                                4107
      if(ix .ne. 0)goto 22321                                              4108
22330 do 22331 j=1,nin                                                     4108
      k=m(j)                                                               4109
      if(xv(k)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 22351                   4109
      ix=1                                                                 4109
      goto 22332                                                           4109
22351 continue                                                             4111
22331 continue                                                             4112
22332 continue                                                             4112
22321 continue                                                             4113
22360 do 22361 i=1,no                                                      4113
      fi=b(0,ic)+g(i,ic)                                                   4115
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)))         4116
      fi=min(max(exmn,fi),exmx)                                            4116
      sxp(i)=sxp(i)-q(i,ic)                                                4117
      q(i,ic)=min(max(emin*sxp(i),exp(fi)),emax*sxp(i))                    4118
      sxp(i)=sxp(i)+q(i,ic)                                                4119
22361 continue                                                             4120
22362 continue                                                             4120
22301 continue                                                             4121
22302 continue                                                             4121
      s=-sum(b(0,:))/nc                                                    4121
      b(0,:)=b(0,:)+s                                                      4122
      if(jx.gt.0)goto 22012                                                4123
      if(ix .ne. 0)goto 22381                                              4124
22390 do 22391 k=1,ni                                                      4124
      if(ixx(k).eq.1)goto 22391                                            4124
      if(ju(k).eq.0)goto 22391                                             4124
      ga(k)=0.0                                                            4124
22391 continue                                                             4125
22392 continue                                                             4125
22400 do 22401 ic=1,nc                                                     4125
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      4126
22410 do 22411 k=1,ni                                                      4126
      if(ixx(k).eq.1)goto 22411                                            4126
      if(ju(k).eq.0)goto 22411                                             4127
      ga(k)=ga(k)+dot_product(r(:,ic),x(:,k))**2                           4128
22411 continue                                                             4129
22412 continue                                                             4129
22401 continue                                                             4130
22402 continue                                                             4130
      ga=sqrt(ga)                                                          4131
22420 do 22421 k=1,ni                                                      4131
      if(ixx(k).eq.1)goto 22421                                            4131
      if(ju(k).eq.0)goto 22421                                             4132
      if(ga(k) .le. al1*vp(k))goto 22441                                   4132
      ixx(k)=1                                                             4132
      ix=1                                                                 4132
22441 continue                                                             4133
22421 continue                                                             4134
22422 continue                                                             4134
      if(ix.eq.1) go to 11020                                              4135
      goto 22012                                                           4136
22381 continue                                                             4137
      goto 22011                                                           4138
22012 continue                                                             4138
      if(kx .le. 0)goto 22461                                              4138
      jerr=-20000-ilm                                                      4138
      goto 21932                                                           4138
22461 continue                                                             4139
      if(jx .le. 0)goto 22481                                              4139
      jerr=-10000-ilm                                                      4139
      goto 21932                                                           4139
22481 continue                                                             4139
      devi=0.0                                                             4140
22490 do 22491 ic=1,nc                                                     4141
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          4141
      a0(ic,ilm)=b(0,ic)                                                   4142
22500 do 22501 i=1,no                                                      4142
      if(y(i,ic).le.0.0)goto 22501                                         4143
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           4144
22501 continue                                                             4145
22502 continue                                                             4145
22491 continue                                                             4146
22492 continue                                                             4146
      kin(ilm)=nin                                                         4146
      alm(ilm)=al                                                          4146
      lmu=ilm                                                              4147
      dev(ilm)=(dev1-devi)/dev0                                            4148
      if(ilm.lt.mnl)goto 21931                                             4148
      if(flmin.ge.1.0)goto 21931                                           4149
      me=0                                                                 4149
22510 do 22511 j=1,nin                                                     4149
      if(a(j,1,ilm).ne.0.0) me=me+1                                        4149
22511 continue                                                             4149
22512 continue                                                             4149
      if(me.gt.ne)goto 21932                                               4150
      if(dev(ilm).gt.devmax)goto 21932                                     4150
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 21932                             4151
21931 continue                                                             4152
21932 continue                                                             4152
      g=log(q)                                                             4152
22520 do 22521 i=1,no                                                      4152
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         4152
22521 continue                                                             4153
22522 continue                                                             4153
      deallocate(sxp,b,bs,r,q,mm,is,ga,ixx,gk,del,sxpl)                    4154
      return                                                               4155
      end                                                                  4156
      subroutine multsprlognetn(parm,no,ni,nc,x,ix,jx,y,g,w,ju,vp,cl,ne,   4158 
     *nx,nlam,  flmin,ulam,shri,intr,maxit,xv,xb,xs,lmu,a0,a,m,kin,dev0,
     *dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   4159
      double precision x(*),y(no,nc),g(no,nc),w(no),vp(ni)                 4160
      double precision ulam(nlam),xb(ni),xs(ni),xv(ni)                     4161
      double precision a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(   4162 
     *2,ni)
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           4163
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
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               4177
      exmn=-exmx                                                           4178
      allocate(mm(1:ni),stat=jerr)                                         4179
      if(jerr.ne.0) return                                                 4180
      allocate(ga(1:ni),stat=jerr)                                         4181
      if(jerr.ne.0) return                                                 4182
      allocate(gk(1:nc),stat=jerr)                                         4183
      if(jerr.ne.0) return                                                 4184
      allocate(del(1:nc),stat=jerr)                                        4185
      if(jerr.ne.0) return                                                 4186
      allocate(iy(1:ni),stat=jerr)                                         4187
      if(jerr.ne.0) return                                                 4188
      allocate(is(1:max(nc,ni)),stat=jerr)                                 4189
      if(jerr.ne.0) return                                                 4190
      allocate(sxp(1:no),stat=jerr)                                        4191
      if(jerr.ne.0) return                                                 4192
      allocate(sxpl(1:no),stat=jerr)                                       4193
      if(jerr.ne.0) return                                                 4194
      allocate(svr(1:nc),stat=jerr)                                        4195
      if(jerr.ne.0) return                                                 4196
      allocate(sc(1:no),stat=jerr)                                         4197
      if(jerr.ne.0) return                                                 4198
      allocate(isc(1:nc),stat=jerr)                                        4199
      if(jerr.ne.0) return                                                 4200
      pmax=1.0-pmin                                                        4200
      emin=pmin/pmax                                                       4200
      emax=1.0/emin                                                        4201
      bta=parm                                                             4201
      omb=1.0-bta                                                          4201
      dev1=0.0                                                             4201
      dev0=0.0                                                             4202
22530 do 22531 ic=1,nc                                                     4202
      q0=dot_product(w,y(:,ic))                                            4203
      if(q0 .gt. pmin)goto 22551                                           4203
      jerr =8000+ic                                                        4203
      return                                                               4203
22551 continue                                                             4204
      if(q0 .lt. pmax)goto 22571                                           4204
      jerr =9000+ic                                                        4204
      return                                                               4204
22571 continue                                                             4205
      b(1:ni,ic)=0.0                                                       4206
      if(intr .ne. 0)goto 22591                                            4206
      q0=1.0/nc                                                            4206
      b(0,ic)=0.0                                                          4206
      goto 22601                                                           4207
22591 continue                                                             4207
      b(0,ic)=log(q0)                                                      4207
      dev1=dev1-q0*b(0,ic)                                                 4207
22601 continue                                                             4208
22581 continue                                                             4208
22531 continue                                                             4209
22532 continue                                                             4209
      if(intr.eq.0) dev1=log(float(nc))                                    4209
      iy=0                                                                 4209
      al=0.0                                                               4210
      if(nonzero(no*nc,g) .ne. 0)goto 22621                                4211
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         4211
      sxp=0.0                                                              4212
22630 do 22631 ic=1,nc                                                     4212
      q(:,ic)=exp(b(0,ic))                                                 4212
      sxp=sxp+q(:,ic)                                                      4212
22631 continue                                                             4213
22632 continue                                                             4213
      goto 22641                                                           4214
22621 continue                                                             4214
22650 do 22651 i=1,no                                                      4214
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         4214
22651 continue                                                             4214
22652 continue                                                             4214
      sxp=0.0                                                              4215
      if(intr .ne. 0)goto 22671                                            4215
      b(0,:)=0.0                                                           4215
      goto 22681                                                           4216
22671 continue                                                             4216
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 4216
      if(jerr.ne.0) return                                                 4216
22681 continue                                                             4217
22661 continue                                                             4217
      dev1=0.0                                                             4218
22690 do 22691 ic=1,nc                                                     4218
      q(:,ic)=b(0,ic)+g(:,ic)                                              4219
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             4220
      q(:,ic)=exp(q(:,ic))                                                 4220
      sxp=sxp+q(:,ic)                                                      4221
22691 continue                                                             4222
22692 continue                                                             4222
      sxpl=w*log(sxp)                                                      4222
22700 do 22701 ic=1,nc                                                     4222
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  4222
22701 continue                                                             4223
22702 continue                                                             4223
22641 continue                                                             4224
22611 continue                                                             4224
22710 do 22711 ic=1,nc                                                     4224
22720 do 22721 i=1,no                                                      4224
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               4224
22721 continue                                                             4224
22722 continue                                                             4224
22711 continue                                                             4225
22712 continue                                                             4225
      dev0=dev0+dev1                                                       4227
      alf=1.0                                                              4229
      if(flmin .ge. 1.0)goto 22741                                         4229
      eqs=max(eps,flmin)                                                   4229
      alf=eqs**(1.0/(nlam-1))                                              4229
22741 continue                                                             4230
      m=0                                                                  4230
      mm=0                                                                 4230
      nin=0                                                                4230
      nlp=0                                                                4230
      mnl=min(mnlam,nlam)                                                  4230
      bs=0.0                                                               4231
      shr=shri*dev0                                                        4231
      ga=0.0                                                               4232
22750 do 22751 ic=1,nc                                                     4232
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      4232
      svr(ic)=sum(r(:,ic))                                                 4233
22760 do 22761 j=1,ni                                                      4233
      if(ju(j).eq.0)goto 22761                                             4234
      jb=ix(j)                                                             4234
      je=ix(j+1)-1                                                         4235
      gj=dot_product(r(jx(jb:je),ic),x(jb:je))                             4236
      ga(j)=ga(j)+((gj-svr(ic)*xb(j))/xs(j))**2                            4237
22761 continue                                                             4238
22762 continue                                                             4238
22751 continue                                                             4239
22752 continue                                                             4239
      ga=sqrt(ga)                                                          4240
22770 do 22771 ilm=1,nlam                                                  4240
      al0=al                                                               4241
      if(flmin .lt. 1.0)goto 22791                                         4241
      al=ulam(ilm)                                                         4241
      goto 22781                                                           4242
22791 if(ilm .le. 2)goto 22801                                             4242
      al=al*alf                                                            4242
      goto 22781                                                           4243
22801 if(ilm .ne. 1)goto 22811                                             4243
      al=big                                                               4243
      goto 22821                                                           4244
22811 continue                                                             4244
      al0=0.0                                                              4245
22830 do 22831 j=1,ni                                                      4245
      if(ju(j).eq.0)goto 22831                                             4245
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            4245
22831 continue                                                             4246
22832 continue                                                             4246
      al0=al0/max(bta,1.0d-3)                                              4246
      al=alf*al0                                                           4247
22821 continue                                                             4248
22781 continue                                                             4248
      al2=al*omb                                                           4248
      al1=al*bta                                                           4248
      tlam=bta*(2.0*al-al0)                                                4249
22840 do 22841 k=1,ni                                                      4249
      if(iy(k).eq.1)goto 22841                                             4249
      if(ju(k).eq.0)goto 22841                                             4250
      if(ga(k).gt.tlam*vp(k)) iy(k)=1                                      4251
22841 continue                                                             4252
22842 continue                                                             4252
11020 continue                                                             4253
22850 continue                                                             4253
22851 continue                                                             4253
      ixx=0                                                                4253
      jxx=ixx                                                              4253
      kxx=jxx                                                              4253
      t=0.0                                                                4254
22860 do 22861 ic=1,nc                                                     4254
      t=max(t,maxval(q(:,ic)*(1.0-q(:,ic)/sxp)/sxp))                       4254
22861 continue                                                             4255
22862 continue                                                             4255
      if(t .ge. eps)goto 22881                                             4255
      kxx=1                                                                4255
      goto 22852                                                           4255
22881 continue                                                             4255
      t=2.0*t                                                              4255
      alt=al1/t                                                            4255
      al2t=al2/t                                                           4256
22890 do 22891 ic=1,nc                                                     4256
      bs(0,ic)=b(0,ic)                                                     4256
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          4257
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)/t                                    4257
      svr(ic)=sum(r(:,ic))                                                 4258
      if(intr .eq. 0)goto 22911                                            4258
      b(0,ic)=b(0,ic)+svr(ic)                                              4258
      r(:,ic)=r(:,ic)-svr(ic)*w                                            4259
      dlx=max(dlx,svr(ic)**2)                                              4260
22911 continue                                                             4261
22891 continue                                                             4262
22892 continue                                                             4262
22920 continue                                                             4262
22921 continue                                                             4262
      nlp=nlp+nc                                                           4262
      dlx=0.0                                                              4263
22930 do 22931 k=1,ni                                                      4263
      if(iy(k).eq.0)goto 22931                                             4264
      jb=ix(k)                                                             4264
      je=ix(k+1)-1                                                         4264
      del=b(k,:)                                                           4264
      gkn=0.0                                                              4265
22940 do 22941 ic=1,nc                                                     4266
      u=(dot_product(r(jx(jb:je),ic),x(jb:je))-svr(ic)*xb(k))/xs(k)        4267
      gk(ic)=u+del(ic)*xv(k)                                               4267
      gkn=gkn+gk(ic)**2                                                    4268
22941 continue                                                             4269
22942 continue                                                             4269
      gkn=sqrt(gkn)                                                        4269
      u=1.0-alt*vp(k)/gkn                                                  4270
      if(u .gt. 0.0)goto 22961                                             4270
      b(k,:)=0.0                                                           4270
      goto 22971                                                           4271
22961 continue                                                             4272
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                     4273
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),cl(2,k),  vp(k)*al2t,alt*vp(   4275 
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                                 4276
22971 continue                                                             4277
22951 continue                                                             4277
      del=b(k,:)-del                                                       4277
      if(maxval(abs(del)).le.0.0)goto 22931                                4278
22980 do 22981 ic=1,nc                                                     4278
      dlx=max(dlx,xv(k)*del(ic)**2)                                        4279
      r(jx(jb:je),ic)=r(jx(jb:je),ic)  -del(ic)*w(jx(jb:je))*(x(jb:je)-x   4281 
     *b(k))/xs(k)
22981 continue                                                             4282
22982 continue                                                             4282
      if(mm(k) .ne. 0)goto 23001                                           4282
      nin=nin+1                                                            4283
      if(nin .le. nx)goto 23021                                            4283
      jxx=1                                                                4283
      goto 22932                                                           4283
23021 continue                                                             4284
      mm(k)=nin                                                            4284
      m(nin)=k                                                             4285
23001 continue                                                             4286
22931 continue                                                             4287
22932 continue                                                             4287
      if(jxx.gt.0)goto 22922                                               4288
      if(dlx.lt.shr)goto 22922                                             4288
      if(nlp .le. maxit)goto 23041                                         4288
      jerr=-ilm                                                            4288
      return                                                               4288
23041 continue                                                             4289
23050 continue                                                             4289
23051 continue                                                             4289
      nlp=nlp+nc                                                           4289
      dlx=0.0                                                              4290
23060 do 23061 l=1,nin                                                     4290
      k=m(l)                                                               4290
      jb=ix(k)                                                             4290
      je=ix(k+1)-1                                                         4290
      del=b(k,:)                                                           4290
      gkn=0.0                                                              4291
23070 do 23071 ic=1,nc                                                     4292
      u=(dot_product(r(jx(jb:je),ic),x(jb:je))  -svr(ic)*xb(k))/xs(k)      4294
      gk(ic)=u+del(ic)*xv(k)                                               4294
      gkn=gkn+gk(ic)**2                                                    4295
23071 continue                                                             4296
23072 continue                                                             4296
      gkn=sqrt(gkn)                                                        4296
      u=1.0-alt*vp(k)/gkn                                                  4297
      if(u .gt. 0.0)goto 23091                                             4297
      b(k,:)=0.0                                                           4297
      goto 23101                                                           4298
23091 continue                                                             4299
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                     4300
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),cl(2,k),  vp(k)*al2t,alt*vp(   4302 
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                                 4303
23101 continue                                                             4304
23081 continue                                                             4304
      del=b(k,:)-del                                                       4304
      if(maxval(abs(del)).le.0.0)goto 23061                                4305
23110 do 23111 ic=1,nc                                                     4305
      dlx=max(dlx,xv(k)*del(ic)**2)                                        4306
      r(jx(jb:je),ic)=r(jx(jb:je),ic)  -del(ic)*w(jx(jb:je))*(x(jb:je)-x   4308 
     *b(k))/xs(k)
23111 continue                                                             4309
23112 continue                                                             4309
23061 continue                                                             4310
23062 continue                                                             4310
      if(dlx.lt.shr)goto 23052                                             4310
      if(nlp .le. maxit)goto 23131                                         4310
      jerr=-ilm                                                            4310
      return                                                               4310
23131 continue                                                             4312
      goto 23051                                                           4313
23052 continue                                                             4313
      goto 22921                                                           4314
22922 continue                                                             4314
      if(jxx.gt.0)goto 22852                                               4315
23140 do 23141 ic=1,nc                                                     4316
      if((b(0,ic)-bs(0,ic))**2.gt.shr) ixx=1                               4317
      if(ixx .ne. 0)goto 23161                                             4318
23170 do 23171 j=1,nin                                                     4318
      k=m(j)                                                               4319
      if(xv(k)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 23191                   4319
      ixx=1                                                                4319
      goto 23172                                                           4319
23191 continue                                                             4321
23171 continue                                                             4322
23172 continue                                                             4322
23161 continue                                                             4323
      sc=b(0,ic)+g(:,ic)                                                   4323
      b0=0.0                                                               4324
23200 do 23201 j=1,nin                                                     4324
      l=m(j)                                                               4324
      jb=ix(l)                                                             4324
      je=ix(l+1)-1                                                         4325
      sc(jx(jb:je))=sc(jx(jb:je))+b(l,ic)*x(jb:je)/xs(l)                   4326
      b0=b0-b(l,ic)*xb(l)/xs(l)                                            4327
23201 continue                                                             4328
23202 continue                                                             4328
      sc=min(max(exmn,sc+b0),exmx)                                         4329
      sxp=sxp-q(:,ic)                                                      4330
      q(:,ic)=min(max(emin*sxp,exp(sc)),emax*sxp)                          4331
      sxp=sxp+q(:,ic)                                                      4332
23141 continue                                                             4333
23142 continue                                                             4333
      s=sum(b(0,:))/nc                                                     4333
      b(0,:)=b(0,:)-s                                                      4334
      if(jxx.gt.0)goto 22852                                               4335
      if(ixx .ne. 0)goto 23221                                             4336
23230 do 23231 j=1,ni                                                      4336
      if(iy(j).eq.1)goto 23231                                             4336
      if(ju(j).eq.0)goto 23231                                             4336
      ga(j)=0.0                                                            4336
23231 continue                                                             4337
23232 continue                                                             4337
23240 do 23241 ic=1,nc                                                     4337
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      4338
23250 do 23251 j=1,ni                                                      4338
      if(iy(j).eq.1)goto 23251                                             4338
      if(ju(j).eq.0)goto 23251                                             4339
      jb=ix(j)                                                             4339
      je=ix(j+1)-1                                                         4340
      gj=dot_product(r(jx(jb:je),ic),x(jb:je))                             4341
      ga(j)=ga(j)+((gj-svr(ic)*xb(j))/xs(j))**2                            4342
23251 continue                                                             4343
23252 continue                                                             4343
23241 continue                                                             4344
23242 continue                                                             4344
      ga=sqrt(ga)                                                          4345
23260 do 23261 k=1,ni                                                      4345
      if(iy(k).eq.1)goto 23261                                             4345
      if(ju(k).eq.0)goto 23261                                             4346
      if(ga(k) .le. al1*vp(k))goto 23281                                   4346
      iy(k)=1                                                              4346
      ixx=1                                                                4346
23281 continue                                                             4347
23261 continue                                                             4348
23262 continue                                                             4348
      if(ixx.eq.1) go to 11020                                             4349
      goto 22852                                                           4350
23221 continue                                                             4351
      goto 22851                                                           4352
22852 continue                                                             4352
      if(kxx .le. 0)goto 23301                                             4352
      jerr=-20000-ilm                                                      4352
      goto 22772                                                           4352
23301 continue                                                             4353
      if(jxx .le. 0)goto 23321                                             4353
      jerr=-10000-ilm                                                      4353
      goto 22772                                                           4353
23321 continue                                                             4353
      devi=0.0                                                             4354
23330 do 23331 ic=1,nc                                                     4355
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          4355
      a0(ic,ilm)=b(0,ic)                                                   4356
23340 do 23341 i=1,no                                                      4356
      if(y(i,ic).le.0.0)goto 23341                                         4357
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           4358
23341 continue                                                             4359
23342 continue                                                             4359
23331 continue                                                             4360
23332 continue                                                             4360
      kin(ilm)=nin                                                         4360
      alm(ilm)=al                                                          4360
      lmu=ilm                                                              4361
      dev(ilm)=(dev1-devi)/dev0                                            4362
      if(ilm.lt.mnl)goto 22771                                             4362
      if(flmin.ge.1.0)goto 22771                                           4363
      me=0                                                                 4363
23350 do 23351 j=1,nin                                                     4363
      if(a(j,1,ilm).ne.0.0) me=me+1                                        4363
23351 continue                                                             4363
23352 continue                                                             4363
      if(me.gt.ne)goto 22772                                               4364
      if(dev(ilm).gt.devmax)goto 22772                                     4364
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 22772                             4365
22771 continue                                                             4366
22772 continue                                                             4366
      g=log(q)                                                             4366
23360 do 23361 i=1,no                                                      4366
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         4366
23361 continue                                                             4367
23362 continue                                                             4367
      deallocate(sxp,b,bs,r,q,mm,is,sc,ga,iy,gk,del,sxpl)                  4368
      return                                                               4369
      end                                                                  4370
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
