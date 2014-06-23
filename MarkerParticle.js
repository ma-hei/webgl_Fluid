
function MarkerParticle(){
    
    this.velocity = new Array(2);
    this.position = new Array(2);
    
    this.velocity[0] = 0;
    this.velocity[1] = 0;
    
    this.position[0] = 0;
    this.position[1] = 0;
    
    this.setPosition = function(x,y){
        
        this.position[0] = x;
        this.position[1] = y;
        
    }
    
    this.setVelocity= function(x,y){
        
        this.velocity[0] = x;
        this.velocity[1] = y;
        
    }
    
    this.getVelocity = function(){
        return this.velocity;
    }
    
    this.getPosition = function(){
        return this.position;
    }
    
}