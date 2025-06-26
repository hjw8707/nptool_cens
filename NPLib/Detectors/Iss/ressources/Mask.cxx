void Mask(){

  unsigned int strips = 128 ;
  double dimension = 127.;
  double active = 125.;
  double pitch = active/strips;
  double width = 0.9;
//
//  double dimension = 100;
//  double active = 97;
//  double pitch = 97/97;
//  double width = 0.5;


  // mm per pixel 
  double scale = 125./12800;
  //pitch in pixel
  unsigned int spitch = pitch/scale;
  unsigned int swidth  = width/scale;
  unsigned int sinter = (pitch - width)/(scale);
  cout << spitch << " " << swidth << " " << sinter << endl;



  // image size
  unsigned int size = dimension/scale;

/* tried by me:
  double dimensionL = 128.42;
  double activeL= 127.9;
  double scaleL = 127.9/12800;
  unsigned int size2 = dimensionL/scaleL;
*/

  cout << "Image size: " << size << "x" << size << endl ;
  double* zargb = new double[size*size];

  TASImage* mask = new TASImage("mask",zargb,size,size,0);
  unsigned int* argb = mask->GetArgbArray();
  unsigned int* argb2 = mask->GetArgbArray();
  unsigned int index = 0;
  double border1 = 0.5*(dimension-active);
  double border2 =(border1+strips*pitch); //(active)/scale+border1;
  unsigned int sborder1=border1/scale;
  unsigned int sborder2=border2/scale; // border2/scale;

  for(unsigned int px = 0 ; px < size ; px++){
    for(unsigned int py = 0 ; py < size ; py++){
      if(px%1000==0)
        cout << "\r" << px << "/" << size << flush; 
      // Compute array index
      index = px * size + py; 
      // Inactive sides
      if(px < sborder1|| py < sborder1 || px > sborder2 || py > sborder2)
        argb[index] = 0xffff0000;
      else{ // strips
        unsigned int coord = px-border1;
        unsigned int nbr = coord/spitch;
       // cout << coord << " " << nbr*spitch+sinter << " " << spitch << " " << sinter << " " << spitch << " " << sinter << endl;
        if(coord<(nbr*spitch+sinter))
        // interstrip
          argb[index] = 0xffff0000+(((nbr-1))<<8)+nbr-2; 
        else if (nbr-1 < strips+1)
          argb[index] = 0xff000000 + nbr-1;
        else
          argb[index] = 0xffff0000;
      }
    }
  }
  mask->WriteImage("mask.png");
  delete[] zargb;
  //mask->Draw();

}
