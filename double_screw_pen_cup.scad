what_to_output = "All"; // [All:All, Screw:Screw only, Left:Left nut only, Right: Right nut only, Top: Top Nut, Bottom: Bottom Nut, vase: Vase Mode Screw, model: Show assembly of parts]
nut_type = "round"; // [round:round, hex:Six Sided]
thread_outer_radius = 35; // [1:100]
thread_inner_radius = 32;  // [1:100]
thread_pitch = 7; // [1:50]
left_starts = 7; // [1:20]
right_starts = 3; // [1:20]
screw_length = 80; //[10:300]
nut_radius = 40; //[10:200]
nut_height = 12; //[5:30]
nut_cut = 3; // [2.5:15]
//Increase this if the nut jams on the screw
nut_thread_tolerance = 0.1; //[0:0.05:1]
//radial subdivision of the thread
thread_resolution = 90; //[10:120]
//Show thread direction on nuts as arrow
show_arrow = true;
//Show thread direction on nuts in text
show_direction = true;
//Show number of starts on nuts
show_starts = true;

$fn=$preview ? 30 : 120;
// REVOLVE2 LIBRARY COMPONENTS INLINED FOR THINGIVERSE CUSTOMIZER
// NOTE: REVOLVE2 IS LICENCED UNDER Creative Commons - Attribution

// Vector generation
function linspace(start,stop, n) = let (step=(stop-start)/n) 
  concat( [ for (i = [0:1:n-0.1]) start+i*step], stop);
function range(start,step,stop) = [ for (i=[start:step:stop]) i ];

// Operations of vectors
function front(v) = v[0];
function back(v) = v[len(v)-1];
function reverse(v) = [for (i=[len(v)-1:-1:0]) v[i]];
function empty(v) = len(v)==0;
function flatten(v) = [ for (a = v) for (b = a) b ] ;
function make_circular(v) = concat(v, [front(v)]);
function pos(v, p, L=-1, R=-1) = //binary_search (v must be sorted)
  L == -1 ? pos(v, p, 0, len(v)-1) :
  L <= R ? let( i = floor((L + R) / 2))
  v[i] < p ? pos(v, p, i + 1, R) :
  v[i] > p ? pos(v, p, L, i - 1) :
  i : -1;
function cum_sum(v, r=[0], i=0) =
  len(r)==len(v) ? r : cum_sum(v, concat(r, back(r)+v[i]), i+1);

// helper functions for dealing with raw_points = [index, x, y, z]
function raw2pt(ipt) = [ipt[1],ipt[2],ipt[3]];
function raw2z(ipt) = ipt[3];
function set_raw_z(ipt,z) = [ipt[0], ipt[1], ipt[2], z];
function raw2indx(ipt) = ipt[0];
function raw_pts2pts (raw_pts) = [ for (col=raw_pts) [for (pt=col) raw2pt(pt) ]];
function raw_pts2indx(raw_pts) = [ for (col=raw_pts) [for (pt=col) raw2indx(pt) ]];

// Geometric operations
function rotate_2D(p, a=0) = concat(
  [cos(a)*p.x-sin(a)*p.y, sin(a)*p.x+cos(a)*p.y],
  [for (i=[2:1:len(p)-1]) p[i]]);
function lookup_2D(pts, z) = 
  let(x=[for (p=pts) [p[2],p[0]]], y=[for (p=pts) [p[2],p[1]]])
    [lookup(z,x), lookup(z,y), z];
function do_scale(r, z, table) = r*max(-1e3, lookup(z,table));

// Manipulations of the list of points
function trim_bottom(raw_pts, val) =
  concat([ for (i=[0:1:len(raw_pts)-2]) if ( raw2z(raw_pts[i+1])>val+1e-9 )
    (raw2z(raw_pts[i])<val) ? 
      concat(raw2indx(raw_pts[i]), lookup_2D([raw2pt(raw_pts[i]),raw2pt(raw_pts[i+1])], val)) : raw_pts[i]
  ], [back(raw_pts)]);
    
function trim_top(raw_pts, val) =
  concat([ for (i=[len(raw_pts)-1:-1:1]) if ( raw2z(raw_pts[i-1])<val-1e-9 ) 
    (raw2z(raw_pts[i])>val) ? 
      concat( raw2indx(raw_pts[i]), lookup_2D([raw2pt(raw_pts[i-1]),raw2pt(raw_pts[i])], val)) : raw_pts[i]
  ], [front(raw_pts)]);

function trim(raw_pts, b, t) = 
  [ for (col=raw_pts) reverse(trim_top(trim_bottom(col, b), t)) ];

function add_bottom_top(raw_pts, minz, maxz) = 
  [for (col=raw_pts) let(b=front(col), t=back(col))
    concat([[b[0]-1,b[1],b[2], minz]], col, [[t[0]+1,t[1],t[2], maxz]]) ];

// Helpers for triangulation
function filter_indx(fcs, indx) = 
  [for (fc=fcs) [for (pt=fc) let(j=pos(indx,pt)) if (j!=-1) j]];

function split_square(v,normal=true) = normal ?
  [[v[0],v[1],v[3]], [v[1],v[2],v[3]]] :
// change the splitting diagonal to reduce face stretching
  [[v[0],v[2],v[3]], [v[0],v[1],v[2]]] ;  

function filter(fcs, indx, normal) = let(tmp=filter_indx(fcs,indx))
  flatten([for (f=tmp) if (len(f)>=3) len(f)==3 ? [f] : split_square(f,normal)]);
    
function make_bottom_top_fcs(indx) = let(lns=[for (i=indx) len(i)])
  [cum_sum(lns,[0]), reverse(cum_sum(lns,[front(lns)-1],1))];
  
function expected_n_fcs(indx, n=0, i=0) = 
  i == len(indx) ? n+2 :
    expected_n_fcs(indx, n+len(indx[i])+(i==0?len(back(indx)):len(indx[i-1]))-2, i+1);

// revolve module
module revolve(profile = [[]], length=0, nthreads=1, scale=1,
               preserve_thread_depth=false, preserve_thread_shape=false,
               $fn=$fn, force_alternate_construction=false) {
                 
  // Extend the profile to accomodate the required length and number of threads
  period = back(profile)[0]-front(profile)[0];
  L = (length == 0) ? period : length;
  profile_ext = let(
    Nz = len(profile), 
    zps = nthreads > 0 ?
      range(-period*(nthreads+1),period,L):
      range(-period,period,L-period*nthreads)
  ) concat(
    [for (zp = zps) for ( np=[0:1:Nz-2]) [ profile[np][0]+zp, profile[np][1]]],
    [[ profile[Nz-1][0]+back(zps), profile[Nz-1][1]]] 
  );

  // Prepare some auxiliary variables
  minZ = min([for (p=profile_ext) p[0]])-1+(nthreads<0 ? period*nthreads : 0);
  maxZ = max([for (p=profile_ext) p[0]])+1+(nthreads>0 ? period*nthreads : 0);
  Nz = len(profile_ext)+1;
  z0 = profile[0][0];
  minR = min([for (p=profile) p[1]]);
  maxR = max([for (p=profile) p[1]]);
  Na = ($fn==undef || $fn<3) ? max(5, maxR*max(abs(scale),1)) : $fn;
  stepa = 360/Na;
  scale_table = [[z0,1],[z0+L, scale]];
  maxS = max(scale,1);
  
  // Compute the vector of points
  raw_pts = add_bottom_top(
    [for (an=[0:1:Na-1]) [for (pn=[1:1:Nz-2])
      let(ai=an*stepa, zi=profile_ext[pn][0]+an*period*nthreads/Na, ri=profile_ext[pn][1],
        r_z = scale==1? [ri,zi] : 
          preserve_thread_shape ?
            rotate_2D([ri,zi], atan(minR*(1-scale)/L)) :
          preserve_thread_depth ?
            [do_scale(minR, zi, scale_table)+ri-minR, zi] :
            [do_scale(ri, zi, scale_table), zi],
        r = r_z[0], z = r_z[1]
      )
      if (r>0) concat(pn+an*Nz, rotate_2D([r, 0, z], ai))]],
    minZ, maxZ);
 
  // Extract the vector of indexes
  raw_indx = raw_pts2indx(raw_pts);

  // Compute the faces
  raw_fcs = concat( 
    [for (col=[0:1:Na-2]) for (row=[0:1:Nz-2]) [
      raw_indx[col][row], raw_indx[col][row+1],
      raw_indx[col+1][row+1], raw_indx[col+1][row] 
    ]],
    nthreads>=0 ?
    concat(
      //last side, 4-point faces
      [for (i=[0:1:Nz-nthreads*(len(profile)-1)-2])
        [raw_indx[Na-1][i], raw_indx[Na-1][i+1], 
         raw_indx[0][i+nthreads*(len(profile)-1)+1],
         raw_indx[0][i+nthreads*(len(profile)-1)]]],
      //bottom connection, 3-point faces
      [for (i=[0:1:nthreads*(len(profile)-1)-1]) 
        [raw_indx[Na-1][0], raw_indx[0][i+1],raw_indx[0][i]]],
      //top connection, 3-point faces
      [for (i=[Nz-nthreads*(len(profile)-1)-1:1:Nz-2]) 
        [raw_indx[Na-1][i], raw_indx[Na-1][i+1],raw_indx[0][Nz-1]]]
      ) :
    concat(
      //last side, 4-point faces
      [for (i=[0:1:Nz-nthreads*(1-len(profile))-2])
        [raw_indx[Na-1][i+nthreads*(1-len(profile))],
         raw_indx[Na-1][i+nthreads*(1-len(profile))+1],
         raw_indx[0][i+1], raw_indx[0][i] ]],
      //bottom connection, 3-point faces
      [for (i=[0:1:nthreads*(1-len(profile))-1]) 
        [raw_indx[Na-1][i], raw_indx[Na-1][i+1],raw_indx[0][0]]],
      //top connection, 3-point faces
      [for (i=[Nz-nthreads*(1-len(profile))-1:1:Nz-2]) 
        [raw_indx[Na-1][Nz-1], raw_indx[0][i+1],raw_indx[0][i]]]
      )
  );

  if (force_alternate_construction) {
    revolve_alternate_construction(raw_pts, raw_fcs, nthreads, maxR*maxS, z0, L);
  } else {
    trim_raw_pts = trim(raw_pts,z0,z0+L);
    trim_indx = raw_pts2indx(trim_raw_pts);
    trim_fcs = concat(filter(raw_fcs, flatten(trim_indx), normal=(nthreads>=0)),
                      make_bottom_top_fcs(trim_indx));

    n_fcs = len(trim_fcs);
    exp_n_fcs = expected_n_fcs(trim_indx);
    if (n_fcs < exp_n_fcs) {
      echo(str("REVOLVE INFO: Using alternate construction."));
      revolve_alternate_construction(raw_pts, raw_fcs, nthreads, maxR*maxS,  z0, L);
    } else {
      trim_pts = flatten(raw_pts2pts(trim_raw_pts));
      polyhedron( points = trim_pts, faces = trim_fcs, convexity=10);
    }
  }
}

module revolve_alternate_construction(raw_pts, raw_fcs, nthreads, maxR, z0, length) {
  indx = raw_pts2indx(raw_pts);
  pts  = flatten(raw_pts2pts(raw_pts));
  fcs  = concat(filter(raw_fcs, flatten(indx), normal=(nthreads>=0)),
               make_bottom_top_fcs(indx));
  intersection() {
    polyhedron( points = pts, faces = fcs, convexity=10);
    translate([0,0,z0]) cylinder(r=maxR+1, h=length);
  } 
}

////////////////////////////////////////// END OF REVOLVE2 LIBRARY
  
// THIS IS THE ACTUAL CODE YOU NEED IF YOU HAVE REVOLVE2 INSTALLED
// THE CODE BELOW IS PUBLIC DOMAIN, NO LICENCE, ENJOY

module cutter(a=5,w=50,h=100) {
    render()
    intersection() {
        rotate([0,0,-a/2]) translate([0,0,-h/2]) cube([w,50,h]);
        rotate([0,0,a/2]) translate([0,-50,-h/2]) cube([w,50,h]);
    }
}

module wrap(radius=nut_radius-1, height=nut_height){
  step=5;
  for( i=[0:step:360+2*step] ) {
    rotate([0,0,i])
    intersection() {
        translate([radius,-i*PI/180*radius,0]) children();
        cutter(step,radius+20,height+20);
    }
  }
}

module direction(threads=1,size=nut_height/2,height=2,upsidedown=false,count=false){
  translate([0,7.5,0])
  rotate([90,0,90])
    linear_extrude(height=height,convexity=2) 
      if ( upsidedown ){
        rotate([180,0,0]) mirror([1,0,0])
        if( count )
          show_count(threads,size);
        else
          direction_internal(threads,size);
      }
      else{
        if( count )
          show_count(threads,size);
        else
          direction_internal(threads,size);
      }
}

module direction_internal(threads=1,size=nut_height/2){
  if( threads < 0 ){
    text("LEFT",size,halign="center",valign="center");
  }
  else{
  text("RIGHT",size,halign="center",valign="center");
  }
}

module show_count(threads=1,size=nut_height/2){
  text(str(abs(threads)),size,halign="center",valign="center");
}

module arrow(threads=1,width=nut_height/2,height=2){
  translate([0,width/3+width/4,0])
  if(threads < 0){
    translate([0,width,width/2]) 
      arrow_int(width=width,height=height);
  }
  else{
    translate([0,width/2,width/2]) 
    rotate([180,0,0]) 
      arrow_int(width=width,height=height);
  }
  translate([0,width/4,width/3])
    rotate([-90,0,0]) 
      arrow_int(width=width/2,height=2);
}

module arrow_single(threads=1,width=nut_height/2,height=2){
  if(threads > 0){
    translate([0,width,width/2]) 
      arrow_int(width=width,height=height);
  }
  else{
    translate([0,width/2,width/2]) 
    rotate([0,0,180]) 
      arrow_int(width=width,height=height);
  }
}

module arrow_int(width=nut_height/2,height=2){
  rotate([90,0,90]) linear_extrude(height=height) union(){
    circle(d=width,$fn=3);
    translate([-width/2,0,0]) 
      square([width,width/3],center=true);
  }
}

module thread(length = 10, n = 1, tol = 0) {
  revolve(
    [[0, thread_inner_radius+tol],
    [thread_pitch/2,thread_outer_radius+tol],
    [thread_pitch, thread_inner_radius+tol]],
    length = length,
    nthreads = n,
    $fn = thread_resolution
  );
}

module screw(head=true) {
  difference(){
  intersection() {
//    union(){
    thread(length = screw_length+nut_height, n = -left_starts);
    thread(length = screw_length+nut_height, n = right_starts);
    }
    // add a bit of conicity to the tip, to help screwing the nuts
//    cylinder( r1 = thread_inner_radius+screw_length+nut_height,
//              r2 = thread_inner_radius*0.9, 
//              h = screw_length+nut_height);
//  }
    translate([0,0,2]) cylinder(r=thread_inner_radius-2,h=screw_length+nut_height);
  }
  // head
  if( head ){
    if( nut_type == "hex" ){
      minkowski() {
        translate([0,0,1]) 
          cylinder(r=nut_radius-2, h=nut_height-2, $fn=6);
        sphere(r=1, $fn=20);
      }
    }
    else if( nut_type == "round" ){
      translate([0,0,nut_height-nut_cut]) 
        cylinder(r1=nut_radius, r2=nut_radius-nut_cut, h=nut_cut);
      translate([0,0,nut_cut]) 
        cylinder(r=nut_radius, h=nut_height-nut_cut*2);
      cylinder(r1=nut_radius-nut_cut, r2=nut_radius, h=nut_cut);
    }
  }
}

module Nut(n_threads, voffset=0, upsidedown=false) {
  difference() {
    if( nut_type == "hex" ) difference(){
      minkowski() {
        translate([0,0,1]) cylinder(r=nut_radius-2, h=nut_height-2, $fn=6);
        sphere(r=1, $fn=20);
      }
      if( show_arrow ) {
        rotate([0,0,30]) 
          translate([(nut_radius*cos(30))-2,0,nut_height/4]) 
            translate([0,-nut_height/3,0]) 
              arrow(threads=n_threads);
      }
      if( show_direction ){
        rotate([0,0,180+30]) 
          translate([(nut_radius*cos(30))-2,-7,nut_height/2]) 
            direction(threads=n_threads, upsidedown=upsidedown);
      }
      if( show_starts ){
        rotate([0,0,150])
          translate([(nut_radius*cos(30))-2,-7,nut_height/2])
            direction(threads=n_threads, upsidedown=upsidedown, count=true);
      }
    }
    else if(nut_type == "round") difference(){
      union(){
        translate([0,0,nut_height-nut_cut]) 
          cylinder(r1=nut_radius, r2=nut_radius-3, h=3);
        translate([0,0,nut_cut]) cylinder(r=nut_radius, h=nut_height-nut_cut*2);
        cylinder(r1=nut_radius-nut_cut, r2=nut_radius, h=nut_cut);
      }
      if( show_arrow ) {
        if( $preview )
          %translate([0,0,nut_cut]) wrap() 
            arrow(threads=n_threads);
        else
          translate([0,0,nut_cut]) wrap() 
            arrow(threads=n_threads);
      }
      if( show_direction ){
        rotate([0,0,180]) translate([0,0,nut_height/2]) 
        if( $preview )
          %wrap() 
            direction(threads=n_threads, upsidedown=upsidedown, size=nut_height/2-1.5);
        else
          wrap() 
            direction(threads=n_threads, upsidedown=upsidedown, size=nut_height/2-1.5);
      }
      if( show_starts ){
        rotate([0,0,90])  translate([0,0,nut_height/2])
        if( $preview )
          %wrap() 
            direction(threads=n_threads, upsidedown=upsidedown, count=true);
        else
          wrap() 
            direction(threads=n_threads, upsidedown=upsidedown, count=true);
      
      }
    }
    translate([0,0,-nut_height/2+voffset]) thread(length = nut_height*2, n = n_threads, tol = nut_thread_tolerance);
  }
}

module top_nut(n_threads){
  
  difference(){
    union(){
      if( nut_type == "hex" )
        cylinder(r=nut_radius-nut_cut, h=nut_cut,$fn=6);
      else
        cylinder(r=nut_radius-nut_cut, h=nut_cut);
      Nut(n_threads,upsidedown=true);
    }
    translate([0,0,-0.01]) cylinder(r1=thread_inner_radius-2+nut_cut, r2=thread_inner_radius-2, h=nut_cut+0.1);
  }
}

module vase(){
  screw(head=false);
}

module right_nut(){
  Nut(right_starts);
}

module left_nut(){
  Nut(-left_starts);
}

module bottom(){
  Nut(right_starts,voffset=nut_height/2+1);
}

if ( what_to_output == "All" || what_to_output == "Screw") screw();
if ( what_to_output == "vase") 
  vase();
if ( what_to_output == "All" || what_to_output == "Right") 
  translate([nut_radius*2+3,0,0]) right_nut();
if ( what_to_output == "All" || what_to_output == "Left") 
  translate([0,-(nut_radius*2+3),0]) left_nut();
if ( what_to_output == "All" || what_to_output == "Top") 
  translate([nut_radius*2+3,-(nut_radius*2+3),0]) top_nut(right_starts);
if ( what_to_output == "Bottom") 
  bottom();
if ( what_to_output == "model" ){
  screw();
  translate([0,0,screw_length/2+nut_height])right_nut();
  translate([0,0,screw_length/2]) left_nut();
  translate([0,0,screw_length+nut_height+nut_cut]) rotate([180,0,0]) top_nut(right_starts);
}