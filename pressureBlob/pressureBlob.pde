ArrayList<terrain> terrainSet = new ArrayList<terrain>();
ArrayList<PVector> tPoints    = new ArrayList<PVector>();
boolean terrainStart = false;
float cDist = 2.0;

ArrayList<blob> blobs = new ArrayList<blob>();
float g = 9.81*0.001;

/*
############
############
. Made check pos+vel INSTEAD OF just pos
- Move to closest P as well as the vel change
############
############
*/

void setup(){
    size(800,800);
    blob newBlob = new blob( new PVector(width/2.0, height/2.0), 32, 100.0 );
    blobs.add(newBlob);
}
void draw(){
    //frameRate(1);
    background(60,60,60);
    drawBlobs();
    displayTerrain();

    for(int i=0; i<blobs.get(0).nodes.size(); i++){
        println(blobs.get(0).nodes.get(i).pos);
    }
    println("---");

    showTpoints();
    overlay();
}
void keyPressed(){
    if(key == '1'){
        blob newBlob = new blob( new PVector(mouseX,mouseY), 16, 30.0 );
        blobs.add(newBlob);
    }
    if(key == '2'){
        blobs.get(0).nodes.get(0).pos.x = mouseX;
        blobs.get(0).nodes.get(0).pos.y = mouseY;
    }
    if(key == '3'){
        //Continue terrain
        if(terrainStart){
            tPoints.add( new PVector(mouseX,mouseY) );
        }
    }
    if(key == '4'){
        //Start / End terrain
        if(!terrainStart){
            tPoints.clear();
            tPoints.add( new PVector(mouseX,mouseY) );
        }
        if(terrainStart){
            terrain newTerrain = new terrain(tPoints);
            terrainSet.add(newTerrain);
            tPoints.clear();
        }
        terrainStart = !terrainStart;
    }
    if(key == 'w'){
      blobs.get(0).qGas+=0.1*8.31*21.0;
    }
    if(key == 's'){
      blobs.get(0).qGas-=0.1*8.31*21.0;
    }
}

void overlay(){
    pushStyle();
    fill(255);
    text(frameRate, 30,30);
    popStyle();
}
void drawBlobs(){
    for(int i=0; i<blobs.size(); i++){
        blobs.get(i).display();
        blobs.get(i).update();
    }
}
void displayTerrain(){
    for(int i=0; i<terrainSet.size(); i++){
        terrainSet.get(i).display();
    }
}
void showTpoints(){
    pushStyle();
    fill(255,0,0);
    stroke(255,255,255);
    for(int i=0; i<tPoints.size(); i++){
        ellipse(tPoints.get(i).x, tPoints.get(i).y, 10, 10);
    }
    popStyle();
}

class terrain{
    ArrayList<face> faces      = new ArrayList<face>();
    ArrayList<PVector> rPoints = new ArrayList<PVector>();
    
    PVector ulBounding = new PVector(0,0); //Upper left corner / lowestX, highestY
    PVector drBounding = new PVector(0,0); //Down right corner / highestX,lowestY

    terrain(ArrayList<PVector> pointCloud){
        copyPointCloud(pointCloud);
        createFaces();
        calcBounding();
    }

    void display(){
        displayFaces();
        displayBounds();
    }
    void displayFaces(){
        pushStyle();
        for(int i=0; i<faces.size(); i++){
            faces.get(i).display();
        }
        popStyle();
    }
    void displayBounds(){
        pushStyle();
        rectMode(CORNERS);
        stroke(237, 223, 66);
        strokeWeight(2);
        noFill();
        rect(ulBounding.x,ulBounding.y, drBounding.x,drBounding.y);
        popStyle();
    }
    void createFaces(){
        /*
        Takes order raw points (rPoints) and converts to a closed polygon of faces, going clockwise
        */
        for(int i=0; i<rPoints.size(); i++){
            face newFace = new face( new PVector(rPoints.get(i).x, rPoints.get(i).y), new PVector(rPoints.get((i+1)%(rPoints.size())).x, rPoints.get((i+1)%(rPoints.size())).y) );
            faces.add(newFace);
        }
    }
    void calcBounding(){
        /*
        Calculates the bounding box for this polygon
        */
        float lbX = rPoints.get(0).x;
        float ubX = rPoints.get(0).x;
        float lbY = rPoints.get(0).y;
        float ubY = rPoints.get(0).y;
        for(int i=0; i<rPoints.size(); i++){
            if(rPoints.get(i).x < lbX){
                lbX = rPoints.get(i).x;
            }
            if(rPoints.get(i).x > ubX){
                ubX = rPoints.get(i).x;
            }
            if(rPoints.get(i).y < lbY){
                lbY = rPoints.get(i).y;
            }
            if(rPoints.get(i).y > ubY){
                ubY = rPoints.get(i).y;
            }
        }
        ulBounding = new PVector(lbX -2.0*cDist, ubY +2.0*cDist);
        drBounding = new PVector(ubX +2.0*cDist, lbY -2.0*cDist);
    }
    void copyPointCloud(ArrayList<PVector> pCloud){
        rPoints.clear();
        for(int i=0; i<pCloud.size(); i++){
            rPoints.add( new PVector(pCloud.get(i).x, pCloud.get(i).y) );
        }
    }
    boolean isWithinBounds(node being){
        /*
        Checks if the entity is within the bounding box area
        ## CAN ADD pos+vel FOR MORE STABILITY ##
        */
        boolean cond1 = (ulBounding.x < being.pos.x) && (being.pos.x < drBounding.x);
        boolean cond2 = (drBounding.y < being.pos.y) && (being.pos.y < ulBounding.y);
        if(cond1 && cond2){
            return true;
        }
        else{
            return false;
        }
    }
    void checkFaceCollision(node being){
        /*
        Checks if an entity is colliding with any of the faces on this piece of terrain
        -If collides
            -Set correct tangent and normal values
            -Reflect entity
        */
        for(int i=0; i<faces.size(); i++){
            float d1 = sqrt( pow(being.pos.x - faces.get(i).pos1.x,2) + pow(being.pos.y - faces.get(i).pos1.y,2) );
            float d2 = sqrt( pow(being.pos.x - faces.get(i).pos2.x,2) + pow(being.pos.y - faces.get(i).pos2.y,2) );
            if(d1+d2 < faces.get(i).l+cDist){    //Colliding with this face
                faces.get(i).recalcNorms(being);
                faces.get(i).findP(being);
                faces.get(i).reboundEntity(being, new PVector(d1,d2), new PVector( sqrt(pow(faces.get(i).pos1.x-faces.get(i).P.x,2)+pow(faces.get(i).pos1.y-faces.get(i).P.y,2)), sqrt(pow(faces.get(i).pos2.x-faces.get(i).P.x,2)+pow(faces.get(i).pos2.y-faces.get(i).P.y,2)) ) );
                break;
            }
        }
    }
}

class face{
    PVector pos1;
    PVector pos2;
    PVector P = new PVector(0,0);  //Closest approach

    PVector tangent = new PVector(0,0);
    PVector normal  = new PVector(0,0);
    float l;
    float e = 0.99;

    face(PVector position1, PVector position2){
        pos1 = position1;
        pos2 = position2;
        l = sqrt( pow(pos1.x-pos2.x,2) + pow(pos1.y-pos2.y,2) );
        tangent = new PVector( (pos2.x-pos1.x) / l, (pos2.y-pos1.y) / l );
        normal  = new PVector(tangent.y, -tangent.x);
    }

    void display(){
        displayFace();
        displayPoints();
        displayNorms();
    }
    void displayFace(){
        //Face itself
        pushStyle();
        stroke(255,255,255);
        strokeWeight(3);
        line(pos1.x,pos1.y, pos2.x,pos2.y);
        popStyle();
    }
    void displayPoints(){
        //End points
        pushStyle();
        fill(200,200,240);
        strokeWeight(1);
        ellipse(pos1.x, pos1.y, 5,5);
        ellipse(pos2.x, pos2.y, 5,5);
        fill(36, 252, 3);
        ellipse(P.x, P.y, 5,5);
        popStyle();
    }
    void displayNorms(){
        //Tangent and normal lines
        pushStyle();
        stroke(255,0,0);
        line( (pos1.x+pos2.x)/2.0, (pos1.y+pos2.y)/2.0, 30.0*tangent.x +(pos1.x+pos2.x)/2.0, 30.0*tangent.y +(pos1.y+pos2.y)/2.0 );
        stroke(0,255,0);
        line( (pos1.x+pos2.x)/2.0, (pos1.y+pos2.y)/2.0, 20.0*normal.x +(pos1.x+pos2.x)/2.0, 20.0*normal.y +(pos1.y+pos2.y)/2.0 );
        popStyle();
    }
    void recalcNorms(node being){
        if(being.vel.x*tangent.x + being.vel.y*tangent.y < 0){
            tangent.x *= -1;
            tangent.y *= -1;
        }
        if(being.vel.x*normal.x + being.vel.y*normal.y > 0){
            //normal.x *= -1;
            //normal.y *= -1;
        }
    }
    void reboundEntity(node being, PVector L, PVector X){
        /*
        l = distance from pos1 to being, and pos2 to being
        x = distance along face to P from pos1, and same for pos2
        */
        PVector cPerp  = new PVector( ((being.vel.x*normal.x)  + (being.vel.y*normal.y))*e*normal.x , ((being.vel.x*normal.x)  + (being.vel.y*normal.y))*e*normal.y );
        PVector cPara  = new PVector( ((being.vel.x*tangent.x) + (being.vel.y*tangent.y))*tangent.x , ((being.vel.x*tangent.x) + (being.vel.y*tangent.y))*tangent.y );
        PVector newVel = new PVector( cPara.x - cPerp.x, cPara.y - cPerp.y );
        float dDist;
        if(pow(L.x,2)+pow(L.y,2) -pow(X.x,2)-pow(X.y,2) > 0){
            dDist = sqrt( (0.5)*( pow(L.x,2)+pow(L.y,2) -pow(X.x,2)-pow(X.y,2) ) );   // PROBLEM WHEN -VE ROOT
        }
        else{
            dDist = 5.0*cDist;
        }
        println(cDist);
        println("DDIST -> ", dDist);
        being.pos = new PVector(P.x +normal.x*dDist, P.y +normal.y*dDist);
        being.vel = new PVector(newVel.x, newVel.y);
    }
    void findP(node being){
        float numerator1 = (being.pos.y - pos1.y)*(tangent.x);
        float numerator2 = (being.pos.x - pos1.x)*(tangent.y);
        float denominator1 = (normal.x)*(tangent.y);
        float denominator2 = (normal.y)*(tangent.x);
        float lambda = (numerator1 - numerator2) / (denominator1 - denominator2);
        P = new PVector(being.pos.x + (lambda)*(normal.x), being.pos.y + (lambda)*(normal.y));
        //dist = sqrt( pow(P.x - pos.x,2) + pow(P.y - pos.y,2) );
    }
}

class node{
    PVector pos;
    PVector vel   = new PVector(0,0);
    PVector force = new PVector(0,0);

    float m;
    float k = 0.1;
    float dConst = 0.01;

    node(PVector initialPos, float mass){
        pos = initialPos;
        m   = mass;
    }

    void display(){
        pushStyle();
        fill(255);
        ellipse(pos.x, pos.y, 10, 10);
        popStyle();
    }
    void calcForce(node n1, node n2, float nLen, float qGas, float volume){
        /*
        n1 and n2 are the nodes cw and ccw of this node
        nLen is the natural length
        Nodes feel;
        .(0)Gravity
        .(1)Elastic force between adjacent
        .(2)Pressure force based on volume
        .(3)Drag forces
        */
        force.x = 0;
        force.y = 0;
        //(0)Gravity
        force.x += 0;
        force.y += m*g;

        //(1)F = kx
        PVector r;
        float rM;
        float x;

        r  = new PVector(n1.pos.x -pos.x, n1.pos.y -pos.y);
        rM = sqrt( pow(r.x,2) + pow(r.y,2) );
        x  = rM - nLen; //+ve means actual in r direction e.g contracting / longer than natural length
        force.x += r.x *(1.0/rM)*(k*x);
        force.y += r.y *(1.0/rM)*(k*x);

        r  = new PVector(n2.pos.x -pos.x, n2.pos.y -pos.y);
        rM = sqrt( pow(r.x,2) + pow(r.y,2) );
        x  = rM - nLen; //+ve means actual in r direction e.g contracting / longer than natural length
        force.x += r.x *(1.0/rM)*(k*x);
        force.y += r.y *(1.0/rM)*(k*x);

        //(2)P = nRT/V, F = PA
        
        PVector tang;
        PVector norm;   //Will all point outwards if one does
        float s1 = sqrt( pow(n1.pos.x -pos.x,2) + pow(n1.pos.y -pos.y,2) );
        tang = new PVector( (n1.pos.x -pos.x)/(s1), (n1.pos.y -pos.y)/(s1) );
        norm = new PVector( -tang.y, tang.x );  //CW
        force.x += (s1)*(qGas / volume) *(norm.x);
        force.y += (s1)*(qGas / volume) *(norm.y);
        float s2 = sqrt( pow(n2.pos.x -pos.x,2) + pow(n2.pos.y -pos.y,2) );
        tang = new PVector( (n2.pos.x -pos.x)/(s2), (n2.pos.y -pos.y)/(s2) );
        norm = new PVector( tang.y, -tang.x );  //CCW => -ve norm
        force.x += (s2)*(qGas / volume) *(norm.x);
        force.y += (s2)*(qGas / volume) *(norm.y);

        //println("Gas F; -> (", (s2)*(qGas / volume) *(norm.x) ,", ", (s2)*(qGas / volume) *(norm.y) ,")");
        
        //(3)Drag
        force.x -= dConst*vel.x;
        force.y -= dConst*vel.y;
    }
    void calcVel(){
        vel.x += force.x / m;
        vel.y += force.y / m;
    }
    void calcPos(){
        pos.x += vel.x;
        pos.y += vel.y;
    }
    void calcCollision(){
        /*
        Checks if the node is colliding with terrain (####OR screen border####)
        */
        //(1)Screen border
        if( (pos.x+vel.x <0) || (pos.x+vel.x >width) ){
            vel = reboundVel(vel, new PVector(0,1));
        }
        if( (pos.y+vel.y <0) || (pos.y+vel.y >height) ){
            vel = reboundVel(vel, new PVector(1,0));
        }

        //(2)Terrain
        for(int i=0; i<terrainSet.size(); i++){
            if( terrainSet.get(i).isWithinBounds(this) ){
                terrainSet.get(i).checkFaceCollision(this);
                break;
            }
        }
    }
    PVector reboundVel(PVector v1, PVector tWall){
        /*
        ####
        ####
        ####
        BASICALLY OBSOLETE -> ONLY FOR BORDER COLLISION
        ####
        ####
        ####
        v1 = approach velocity
        tWall = wall tangent vector
        v2 = exit velocity
        */
        float e = 0.8;
        PVector nWall = new PVector(tWall.y, -tWall.x);
        if( v1.x*nWall.x + v1.y*nWall.y > 0 ){         //Corrects normal direction
            nWall.x *= -1;nWall.y *= -1;}
        if( v1.x*tWall.x + v1.y*tWall.y <= 0 ){       //Corrects tangent direction
            tWall.x *= -1;tWall.y *= -1;}
        PVector cPerp = new PVector( ((v1.x*nWall.x) + (v1.y*nWall.y))*e*nWall.x , ((v1.x*nWall.x) + (v1.y*nWall.y))*e*nWall.y );
        PVector cPara = new PVector( ((v1.x*tWall.x) + (v1.y*tWall.y))  *tWall.x , ((v1.x*tWall.x) + (v1.y*tWall.y))  *tWall.y );
        PVector v2    = new PVector( cPara.x - cPerp.x, cPara.y - cPerp.y );
        return v2;
    }
}

class blob{
    ArrayList<node> nodes = new ArrayList<node>();

    PVector cPos;

    int n;

    float r;
    float dTheta;
    float l;    //Natural length
    float qGas = 0.05*8.31*21.0;  //Quantity of gas, (nRT representation)

    blob(PVector centralInitialPosition, int nodeNumber, float radius){
        cPos = centralInitialPosition;
        n = nodeNumber;
        r = radius;
        dTheta = (2.0*PI)/(n);
        l = sqrt( (2.0*pow(r,2))*(1.0-cos(dTheta)) );
        createNodes();
    }

    void display(){
        displayBorder();
        //displayNodes();
    }
    void update(){
        updateNodeVals();
    }
    void displayBorder(){
        pushStyle();
        noFill();
        //stroke(149, 232, 225);
        stroke(255*(nodes.get(0).pos.x/width), 100 , 255*(nodes.get(0).pos.y/height));
        strokeWeight(3);
        for(int i=0; i<nodes.size(); i++){
            curve(nodes.get(i %(nodes.size())).pos.x,nodes.get(i %(nodes.size())).pos.y,        nodes.get((i+1)%(nodes.size())).pos.x,nodes.get((i+1)%(nodes.size())).pos.y,  
                  nodes.get((i+2)%(nodes.size())).pos.x,nodes.get((i+2)%(nodes.size())).pos.y,  nodes.get((i+3)%(nodes.size())).pos.x,nodes.get((i+3)%(nodes.size())).pos.y);
        }
        popStyle();
    }
    void displayNodes(){
        for(int i=0; i<nodes.size(); i++){
            nodes.get(i).display();
        }
    }
    void updateNodeVals(){
        for(int i=0; i<nodes.size(); i++){
            if(i != 0){
                nodes.get(i).calcForce( nodes.get((i-1)%(nodes.size())), nodes.get((i+1)%(nodes.size())), l, qGas, findVolume() );
            }
            else{
                nodes.get(i).calcForce( nodes.get(nodes.size()-1), nodes.get(1), l, qGas, findVolume() );
            }
            nodes.get(i).calcVel();
            nodes.get(i).calcCollision();
        }
        for(int i=0; i<nodes.size(); i++){
            nodes.get(i).calcPos();
        }
    }
    void createNodes(){
        /*
        Creates a node ring
        */
        nodes.clear();
        for(int i=0; i<n; i++){
            node newNode = new node( new PVector(cPos.x+ r*cos(i*dTheta), cPos.y+ r*sin(i*dTheta)), 1.0 );
            nodes.add(newNode);
        }
    }
    float findVolume(){
        /*
            X   Y
        A   1   4   e.g (1*9 + 3*2 + 2*3) - (4*3 + 9*2 + 2*7) = 2*Area
        B   3   9
        C   2   2
        D   7   3
        ...
        */
        float volume = 0.0;
        /*
        //SHOELACE THEOREM
        for(int i=0; i<nodes.size()-1; i++){
            volume += (nodes.get(i).pos.x)*(nodes.get(i+1).pos.y) - (nodes.get(i).pos.y)*(nodes.get(i+1).pos.x);
        }
        println("VOLUME -> ",abs(volume/2.0));
        return abs(volume/2.0);
        */
        //CENTRE SUMMED AREAS, A = absinC/2
        recalcCPos();
        for(int i=0; i<nodes.size(); i++){
            PVector r0 = new PVector( nodes.get(i).pos.x -cPos.x,                    nodes.get(i).pos.y -cPos.y );
            PVector r1 = new PVector( nodes.get((i+1)%(nodes.size())).pos.x -cPos.x, nodes.get((i+1)%(nodes.size())).pos.y -cPos.y );
            float mR0 = sqrt( pow(r0.x,2) + pow(r0.y,2) );
            float mR1 = sqrt( pow(r1.x,2) + pow(r1.y,2) );
            float cosTheta = (r0.x*r1.x + r0.y*r1.y) / (mR0*mR1);
            float sinTheta = sqrt( 1 - pow(cosTheta,2) );
            volume += (0.5)*(mR0*mR1)*(sinTheta);
        }
        return volume;
    }
    void recalcCPos(){
        cPos = new PVector(0,0);
        for(int i=0; i<nodes.size(); i++){
            cPos.x += nodes.get(i).pos.x;
            cPos.y += nodes.get(i).pos.y;
        }
        cPos.x /= nodes.size();
        cPos.y /= nodes.size();
    }
}
class yBlob extends blob{
    //pass

    yBlob(PVector cPos, int n, float r){
        super(cPos, n, r);
        //pass
    }

    void display(){
        //pass
    }
}
