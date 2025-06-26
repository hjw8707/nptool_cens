void MaskBack(){
  
  unsigned int strips = 32 ;
  double dimension = 67.975;
  double active = 63.96;
  double pitch = 2;
  double width = 1.960;
  
  // mm per pixel 
  double scale = 0.005;
  //pitch in pixel
  unsigned int spitch = pitch/scale;
  unsigned int swidth  = width/scale;
  unsigned int sinter = spitch - swidth;
  cout << spitch << " " << swidth << " " << sinter << endl;
  
  // image size
  unsigned int size = dimension/scale;
  cout << "Image size: " << size << "x" << size << endl;
  double* zargb = new double[size*size];
  TASImage* mask = new TASImage("mask",zargb,size,size,0);
  unsigned int* argb = mask->GetArgbArray();
  unsigned int index = 0;
  double border1 = 0.5*(dimension-active);
  double border2 = (border1+active); //(active)/scale+border1;
  unsigned int sborder1=402; //border1/scale;
  unsigned int sborder2=size-401; //border2/scale;
  cout << "Border 1 : " << sborder1 << " Border 2 : " << sborder2 << endl;
  
  for(unsigned int py = 0 ; py < size ; py++){
    bool test = true;
    for(unsigned int px = 0 ; px < size ; px++){
      if(px%1000==0) cout << "\r" << px << "/" << size << flush; 
      // Compute array index
      index = px * size + py; 
      // Inactive sides
      if(px < sborder1|| py < sborder1 || px > sborder2 || py > sborder2) argb[index] = 0xffff0000;
      else{ // strips
        unsigned int coord = py-sborder1;
        unsigned int nbr = coord/spitch;
        if(coord>=((nbr+1)*spitch-sinter)){
          // interstrip
          argb[index] = 0xffff0000+((nbr+2)<<8)+nbr+1;          
          // if(test) {
          //   // cout << coord << " " << nbr << " " << nbr*spitch+sinter << " " << spitch << " " << sinter << " " << sborder1 << " " << px << endl;
          //   cout << coord << " " << nbr << " " << hex << argb[index] << dec << endl;
          //   test = false;
          // }
        }
        else if (nbr < strips+1) argb[index] = 0xff000000 + nbr + 1;
        else argb[index] = 0xff00ff00;
      }
    }
  }
  
  
  mask->WriteImage("maskBack.png");
  delete[] zargb;
  mask->Draw();
  
}
