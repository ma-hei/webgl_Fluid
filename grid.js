
function grid(gridSize){
    this.gridSize = gridSize;
    this.timestep = 0.8;
    this.displaySize = 80;
    this.uvalues = new Array((gridSize+1)*gridSize);
    this.vvalues = new Array((gridSize+1)*gridSize);
    this.advecteduValues = new Array((gridSize+1)*gridSize);
    this.advectedvValues = new Array((gridSize+1)*gridSize);
    this.upositionsx = new Array((gridSize+1)*gridSize);
    this.upositionsy = new Array((gridSize+1)*gridSize);
    this.vpositionsx = new Array((gridSize+1)*gridSize);
    this.vpositionsy = new Array((gridSize+1)*gridSize);
    this.pressure = new Array(this.gridSize);
    this.markers = new Array(this.gridSize*this.gridSize*4);
    this.oldVelocityVerticalTemp = new Array(2);
    this.oldVelocityHorizontalTemp = new Array(2);
    
    this.oldPositionVerticalTemp = new Array(2);
    this.oldPositionHorizontalTemp = new Array(2);
    
    this.velocityAtPosition = new Array(2);
    
    this.markerPositionTemp = new Array(2);
    this.rhs = new Array(this.gridSize*this.gridSize);
    
    this.m = new sparse_matrix(this.gridSize*this.gridSize, 5);
    this.veloUpLeft = new Array(2);
    this.veloUpRight = new Array(2);
    this.veloDownLeft = new Array(2);
    this.veloDownRight = new Array(2);
    

    
    for (var i=0; i<this.gridSize*this.gridSize*4;i++){
        this.markers[i] = new MarkerParticle();
    }
    
    for (var i=0;i<this.gridSize;i++){
        this.pressure[i] = new Array(this.gridSize);
    }
    
    for (var i=0;i<this.gridSize;i++){
        for (var k=0;k<this.gridSize;k++){
            this.pressure[i][k] = 0;
        }
    }
    
    for (var i=0; i<(gridSize+1)*gridSize;i++){
        this.uvalues[i]=0;
        this.vvalues[i]=0;
        this.advecteduValues[i]=0;
        this.advectedvValues[i]=0;
    }
    
    this.uvalues[10] = 15;
    
    var cellWidth = this.displaySize / gridSize;
    var tempx = 0;
    var tempy = 0.5*cellWidth;
    for (var i = 0; i<gridSize; i++){
        tempx = 0;
        for (var k=0; k<gridSize+1;k++){
            
            this.upositionsx[i*(gridSize+1)+k] = tempx;
            this.upositionsy[i*(gridSize+1)+k] = tempy;
            this.vpositionsx[i*(gridSize+1)+k] = tempy;
            this.vpositionsy[i*(gridSize+1)+k] = tempx;
            tempx += cellWidth;
        }
        tempy+= cellWidth;
    }
    
    tempx = 0;
    tempy = 0;
    var counter = 0;
    for (var i =0; i<this.gridSize; i++){
        tempx = 0;
        for (var k=0; k<this.gridSize; k++){
            this.markers[counter].setPosition(tempx+0.5*cellWidth-0.25*cellWidth,tempy+0.5*cellWidth-0.25*cellWidth);
            this.markers[counter].setVelocity(0,0);
            counter++;
            this.markers[counter].setPosition(tempx+0.5*cellWidth+0.25*cellWidth,tempy+0.5*cellWidth-0.25*cellWidth);
            this.markers[counter].setVelocity(0,0);
            counter++;
            this.markers[counter].setPosition(tempx+0.5*cellWidth-0.25*cellWidth,tempy+0.5*cellWidth+0.25*cellWidth);
            this.markers[counter].setVelocity(0,0);
            counter++;
            this.markers[counter].setPosition((tempx+0.5*cellWidth+0.25*cellWidth),(tempy+0.5*cellWidth+0.25*cellWidth));
            this.markers[counter].setVelocity(0,0);
            //document.write(counter+"</br>");
            //document.write((tempx+0.5*cellWidth+0.25*cellWidth)+"--"+(tempy+0.5*cellWidth+0.25*cellWidth)+"</br>");
            counter++;
            tempx+=cellWidth;
            
        }
        tempy+=cellWidth;
    }
 
    
    
    this.getVelocityAtVerticalCellFace = function(i, k, oldVelocity){
        
        var velx;
        var vely = 0;
        
        velx = this.uvalues[i*(this.gridSize+1)+k];
        
        if (k==0){
            vely += this.vvalues[i];
            vely += this.vvalues[i+1];
            vely = vely/4;
        }
        else if (k == this.gridSize){
            vely += this.vvalues[(k-1) * (this.gridSize + 1) + i];
            vely += this.vvalues[(k-1) * (this.gridSize + 1) + i + 1];
            vely = vely/4;
        }
        else{
            vely += this.vvalues[(k-1)*(this.gridSize + 1)+i];
            vely += this.vvalues[(k-1)*(this.gridSize + 1) + i + 1];
            vely += this.vvalues[k*(this.gridSize + 1) + i];
            vely += this.vvalues[k*(this.gridSize + 1) + i +1];
            vely = vely/4;
        }
        
        oldVelocity[0] = velx;
        oldVelocity[1] = vely;
        this.oldVelocityVerticalTemp[0] = velx;
        this.oldVelocityVerticalTemp[1] = vely;
        
    }
    
    
    
    this.getOldPositionOfParticleAtVerticalCellFace = function(i,k){
        
        var oldVelocity = new Array(2);
        this.getVelocityAtVerticalCellFace(i,k,oldVelocity);
       
        //oldPosition[0] = this.upositionsx[i*(this.gridSize+1)+k]-this.oldVelocityVerticalTemp[0]*this.timestep;
        //oldPosition[1] = this.upositionsy[i*(this.gridSize+1)+k]-this.oldVelocityVerticalTemp[1]*this.timestep;
        
        this.oldPositionVerticalTemp[0] = this.upositionsx[i*(this.gridSize+1)+k]-this.oldVelocityVerticalTemp[0]*this.timestep;
        this.oldPositionVerticalTemp[1] = this.upositionsy[i*(this.gridSize+1)+k]-this.oldVelocityVerticalTemp[1]*this.timestep;
        oldVelocity = null;
        
        
    }
    
    this.getVelocityAtHorizontalCellFace = function(i,k,oldVelocity){
        
        var vely = this.vvalues[i*(this.gridSize+1)+k];
        var velx = 0;
        
        if (k==0){
            velx = this.uvalues[k*(this.gridSize+1)+i];
            velx += this.uvalues[k*(this.gridSize+1)+i+1];
            velx = velx/4;
        }
        else if (k==this.gridSize){
            velx = this.uvalues[(k-1)*(this.gridSize+1)+i];
            velx += this.uvalues[(k-1)*(this.gridSize+1)+i+1];
            velx = velx/4;
        }
        else{
            velx += this.uvalues[(k-1)*(this.gridSize+1)+i];
            velx += this.uvalues[(k-1)*(this.gridSize+1)+i+1];
            velx += this.uvalues[k*(this.gridSize+1)+i];
            velx += this.uvalues[k*(this.gridSize+1)+i+1];
            velx = velx/4;
            
        }
        
        oldVelocity[0] = velx;
        oldVelocity[1] = vely;
        this.oldVelocityHorizontalTemp[0] = velx;
        this.oldVelocityHorizontalTemp[1] = vely;

        
    }
    
    this.getOldPositionOfParticleAtHorizontalCellFace = function(i,k){
        
        var oldVelocity = new Array(2);
        this.getVelocityAtHorizontalCellFace(i,k,oldVelocity);
        
        //oldPosition[0] = this.vpositionsx[i*(this.gridSize+1)+k] - this.oldVelocityHorizontalTemp[0] * this.timestep;
        //oldPosition[1] = this.vpositionsy[i*(this.gridSize+1)+k] - this.oldVelocityHorizontalTemp[1] * this.timestep;
        
        this.oldPositionHorizontalTemp[0] = this.vpositionsx[i*(this.gridSize+1)+k] - this.oldVelocityHorizontalTemp[0] * this.timestep;
        this.oldPositionHorizontalTemp[1] = this.vpositionsy[i*(this.gridSize+1)+k] - this.oldVelocityHorizontalTemp[1] * this.timestep;

        
        oldVelocity = null;
    }
    
    this.getCellOfPosition = function(x,y,cell){
        
        var cellWidth = this.displaySize / this.gridSize;
        var faceleft = Math.floor(x/cellWidth);
        var facelower = Math.floor(y/cellWidth);
        
        if (faceleft == this.gridSize){
            faceleft=faceleft-1;
        }
        if (facelower ==this.gridSize){
            facelower = facelower-1;
        }
        
        cell[0]=faceleft;
        cell[1]=facelower;
        
        
    }
    
    this.getSectorOfPosition = function(leftFacex, lowerFacey, x, y){
        
        var sector;
        var cellWidth = this.displaySize/this.gridSize;
        var centerx = leftFacex + 0.5*cellWidth;
        var centery = lowerFacey + 0.5*cellWidth;
        
        if (x<=centerx && y<=centery){
            sector = 1;
        }
        else if(x>centerx && y<=centery){
            sector = 2;
        }
        else if(x>centerx && y>centery){
            sector = 3;
        }
        else if(x<=centerx && y>centery){
            sector = 4;
        }
        return sector;
        
    }
    
    this.getVelocityAtCellPosition = function(x,y,which,velo){
        
        velo[0] = -1;
        velo[1] = -1;
        
        if (which == 1){
            if (x==0){
                if (y==0){
                    
                    velo[0] = 0.5*this.uvalues[0];
                    velo[1] = 0.5*this.vvalues[0];
                }
                else{
                    velo[0] = 0.5*(this.uvalues[(y-1)*(this.gridSize)+1]+this.uvalues[y*(this.gridSize+1)]);
                    velo[1] = 0.5*(this.vvalues[y]);
                }
            }
            else{
                if (y==0){
                    velo[0] = 0.5*(this.uvalues[x]);
                    velo[1] = 0.5*(this.vvalues[x*(this.gridSize+1)] + this.vvalues[(x-1)*(this.gridSize+1)]);
                }
                else{
                    velo[0] = 0.5*(this.uvalues[(y-1)*(this.gridSize+1)+x] + this.uvalues[y*(this.gridSize+1)+x]);
                    velo[1] = 0.5*(this.vvalues[x*(this.gridSize+1)+y]+this.vvalues[(x-1)*(this.gridSize+1)+y]);
                }
            }
        }
        else if (which==2){
            if (x==(this.gridSize-1)){
                if (y==0){
                    velo[0] = 0.5*this.uvalues[this.gridSize];
                    velo[1] = 0.5*this.vvalues[(this.gridSize-1)*(this.gridSize+1)];
                }
                else{
                    velo[0] = 0.5*(this.uvalues[y*(this.gridSize+1)+this.gridSize]+this.uvalues[(y-1)*(this.gridSize+1)+this.gridSize]);
                    velo[1] = 0.5*this.vvalues[(this.gridSize-1)*(this.gridSize+1)+y];
                }
            }
            else{
                if (y==0){
                    velo[0] = 0.5*this.uvalues[x+1];
                    velo[1] = 0.5*(this.vvalues[x*(this.gridSize+1)]+this.vvalues[(x+1)*(this.gridSize+1)]);
                }
                else{
                    velo[0] = 0.5*(this.uvalues[(y-1)*(this.gridSize+1)+x+1]+this.uvalues[y*(this.gridSize+1)+x+1]);
                    velo[1] = 0.5*(this.vvalues[x*(this.gridSize+1)+y]+this.vvalues[(x+1)*(this.gridSize+1)+y]);
                }
            }
        }
        else if (which==3){
            if (x==(this.gridSize-1)){
                if (y==this.gridSize-1){
                    velo[0]=0.5*(this.uvalues[(this.gridSize-1)*(this.gridSize+1)+this.gridSize]);
                    velo[1]=0.5*(this.vvalues[(this.gridSize-1)*(this.gridSize+1)+this.gridSize]);
                }
                else{
                    velo[0]=0.5*(this.uvalues[y*(this.gridSize+1)+this.gridSize]+this.uvalues[(y+1)*(this.gridSize+1)+this.gridSize]);
                    velo[1]=0.5*(this.vvalues[(this.gridSize-1)*(this.gridSize+1)+y+1]);
                }
                
            }
            else{
                if (y==(this.gridSize-1)){
                    velo[0]=0.5*(this.uvalues[(this.gridSize-1)*(this.gridSize+1)+x+1]);
                    velo[1]=0.5*(this.vvalues[x*(this.gridSize+1)+this.gridSize]+this.vvalues[(x+1)*(this.gridSize+1)+this.gridSize]);
                }
                else{
                    velo[0]=0.5*(this.uvalues[y*(this.gridSize+1)+x+1]+this.uvalues[(y+1)*(this.gridSize+1)+x+1]);
                    velo[1]=0.5*(this.vvalues[x*(this.gridSize+1)+y+1]+this.vvalues[(x+1)*(this.gridSize+1)+y+1]);
                }
                
            }
        }
        else if (which==4){
            if (x==0){
                if (y==(this.gridSize-1)){
                    velo[0]=0.5*(this.uvalues[(this.gridSize-1)*(this.gridSize+1)]);
                    velo[1]=0.5*(this.vvalues[this.gridSize]);
                }
                else{
                    velo[0]=0.5*(this.uvalues[y*(this.gridSize+1)]+this.uvalues[(y+1)*(this.gridSize+1)]);
                    velo[1]=0.5*(this.vvalues[y+1]);
                }
            }
            else{
                if (y==(this.gridSize-1)){
                    velo[0]=0.5*(this.uvalues[(this.gridSize-1)*(this.gridSize+1)+x]);
                    velo[1]=0.5*(this.vvalues[(x-1)*(this.gridSize+1)+this.gridSize]+this.vvalues[x*(this.gridSize+1)+this.gridSize]);
                }
                else{
                    velo[0]=0.5*(this.uvalues[y*(this.gridSize+1)+x+1]+this.uvalues[(y+1)*(this.gridSize+1)+x+1]);
                    velo[1]=0.5*(this.vvalues[x*(this.gridSize+1)+y+1]+this.vvalues[(x-1)*(this.gridSize+1)+y+1]);
                }
            }
        }
        else if (which==5){
            velo[0]=0.5*(this.uvalues[y*(this.gridSize+1)+x]+this.uvalues[y*(this.gridSize+1)+x+1]);
            velo[1]=0.5*(this.vvalues[x*(this.gridSize+1)+y]+this.vvalues[x*(this.gridSize+1)+y+1]);
        }
        
    }
    
    this.getVelocityAtPosition2 = function(x,y){
        
        
        var cell = new Array(2);
        
        this.getCellOfPosition(x,y,cell);
        //document.write("cell: "+cell[0]+"; "+cell[1]+"</br>");
        
        var sector = this.getSectorOfPosition(this.upositionsx[cell[1]*(this.gridSize+1)+cell[0]], this.vpositionsy[cell[0]*(this.gridSize+1)+cell[1]],x,y);
       // document.write("sector: "+sector+"</br>");
        if (sector==1){
            this.getVelocityAtVerticalCellFace(cell[1],cell[0],this.veloUpLeft);
            this.getVelocityAtCellPosition(cell[0],cell[1],5,this.veloUpRight);
            this.getVelocityAtCellPosition(cell[0],cell[1],1,this.veloDownLeft);
            this.getVelocityAtHorizontalCellFace(cell[0],cell[1],this.veloDownRight);
           
        }
        else if (sector==2){
            this.getVelocityAtCellPosition(cell[0],cell[1],5,this.veloUpLeft);
            this.getVelocityAtVerticalCellFace(cell[1],cell[0]+1,this.veloUpRight);
            this.getVelocityAtHorizontalCellFace(cell[0],cell[1],this.veloDownLeft);
            this.getVelocityAtCellPosition(cell[0],cell[1],2,this.veloDownRight);
        }
        else if (sector==3){
            this.getVelocityAtHorizontalCellFace(cell[0],cell[1]+1,this.veloUpLeft);
            this.getVelocityAtCellPosition(cell[0],cell[1],3,this.veloUpRight);
            this.getVelocityAtCellPosition(cell[0],cell[1],5,this.veloDownLeft);
            this.getVelocityAtVerticalCellFace(cell[1],cell[0]+1,this.veloDownRight);
        }
        else if (sector==4){
            this.getVelocityAtCellPosition(cell[0],cell[1],4,this.veloUpLeft);
            this.getVelocityAtHorizontalCellFace(cell[0],cell[1]+1, this.veloUpRight);
            this.getVelocityAtVerticalCellFace(cell[1],cell[0],this.veloDownLeft);
            this.getVelocityAtCellPosition(cell[0],cell[1],5,this.veloDownRight);
        }
        
        
        this.velocityAtPosition[0] = 0.25*(this.veloUpLeft[0]+this.veloUpRight[0]+this.veloDownLeft[0]+this.veloDownRight[0]);
        this.velocityAtPosition[1] = 0.25*(this.veloUpLeft[1]+this.veloUpRight[1]+this.veloDownLeft[1]+this.veloDownRight[1]);

        
                
        
    }
    
    
    
    this.advectAll = function(){
        
        
        for (var i=0;i<this.gridSize;i++){
            for (var k=0;k<this.gridSize+1;k++){
                this.getOldPositionOfParticleAtVerticalCellFace(i,k);
                this.getVelocityAtPosition2(this.oldPositionVerticalTemp[0], this.oldPositionVerticalTemp[1]);
                this.advecteduValues[i*(this.gridSize+1)+k] = this.velocityAtPosition[0];
                
                this.getOldPositionOfParticleAtHorizontalCellFace(i,k);
                this.getVelocityAtPosition2(this.oldPositionHorizontalTemp[0], this.oldPositionHorizontalTemp[1]);
                this.advectedvValues[i*(this.gridSize+1)+k] = this.velocityAtPosition[1];
                
            }
        }
        
    }
    
    var set = 0;
    
    this.solvePressure = function(){
        
        //m = new sparse_matrix(this.gridSize*this.gridSize, 5);
        
        var hauptKoeff = 4;
        var temp;
        var loesung;
        var indice;
        
        for (var i = 0; i< this.gridSize; i++){
            for (var k = 0; k< this.gridSize; k++){
                
                hauptKoeff = 4;
                loesung = 0;
                indice = i*this.gridSize+k;
                
                if ((k-1) < 0){
                    hauptKoeff-=1;
                }
                else{
                    temp = indice -1;
                    if (set == 0){
                        this.m.set_element(indice, temp, -1);
                    }
                    loesung += this.advecteduValues[i*(this.gridSize+1)+k];
                }
                
                if ((i-1) <0){
                    hauptKoeff-=1;
                }
                else{
                    temp = (i-1)*this.gridSize+k;
                    if (set == 0){
                        this.m.set_element(indice, temp, -1);
                    }
                    loesung += this.advectedvValues[k*(this.gridSize+1)+i];
                }
                
                if ((k+1) > (this.gridSize -1)){
                    hauptKoeff-=1;
                }
                else{
                    temp = indice +1;
                    if (set == 0){
                        this.m.set_element(indice, temp, -1);
                    }
                    loesung -= this.advecteduValues[i*(this.gridSize+1)+k+1];
                }
                
                if ((i+1)>(this.gridSize-1)){
                    hauptKoeff-=1;
                }
                else{
                    temp = (i+1)*this.gridSize + k;
                    if (set == 0){
                        this.m.set_element(indice, temp, -1);
                    }
                    loesung -= this.advectedvValues[k*(this.gridSize+1)+i+1];
                }
                
                temp = indice;
                if (set == 0){
                    this.m.set_element(indice, indice, hauptKoeff);
                }
                this.rhs[indice] = loesung;
            }
        }
        
      
        z = conjugate_gradient2(this.m,this.rhs,40,this.gridSize);
       
        var counter = 0;
        var counter2;
        for (var i=0; i<this.gridSize; i++){
            counter2 = i*this.gridSize;
            for (var k=0; k<this.gridSize;k++){
               
                counter = counter2 + k;
                this.pressure[i][k] = z[counter];
               
                
            }
        }
        
        z = null;
        
        set = 1;
        
    }
    
    this.updateVelocities = function(){
        
        for (var i=0;i<this.gridSize;i++){
            for (var k=1;k<this.gridSize;k++){
                this.uvalues[i*(this.gridSize+1)+k] = this.advecteduValues[i*(this.gridSize+1)+k]-this.pressure[i][k]+this.pressure[i][k-1];
            }
        }
        
        for (var i=0;i<this.gridSize;i++){
            for (var k=1;k<this.gridSize;k++){
                this.vvalues[i*(this.gridSize+1)+k] = this.advectedvValues[i*(this.gridSize+1)+k]-this.pressure[k][i]+this.pressure[k-1][i];
            }
        }
     
        
        
    }
    
    
    
    this.moveMarkers = function(){
        
        
        
        for (var i=0; i<this.gridSize*this.gridSize*4;i++){
            
            this.getVelocityAtPosition2(this.markers[i].getPosition()[0], this.markers[i].getPosition()[1]);
            this.markers[i].setVelocity(this.velocityAtPosition[0], this.velocityAtPosition[1]);
            this.markerPositionTemp[0] = this.markers[i].getPosition()[0] + this.timestep*this.velocityAtPosition[0];
            this.markerPositionTemp[1] = this.markers[i].getPosition()[1] + this.timestep*this.velocityAtPosition[1];
            
            if (this.markerPositionTemp[0] < 0){
                this.markerPositionTemp[0] = 0;
            }
            if (this.markerPositionTemp[1] < 0){
                this.markerPositionTemp[1] = 0;
            }
            
            this.markers[i].setPosition(this.markerPositionTemp[0],this.markerPositionTemp[1]);
            
        }
        
    }
    
    
}