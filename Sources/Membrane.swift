import Foundation
import MathTools

// this handles special case of a membrane

struct Edge : Hashable {
    var v1:Int
    var v2:Int

    func hash(into hasher: inout Hasher) {
        hasher.combine(v1)
        hasher.combine(v2)
    }
}

extension Range {
    func contains(otherRange: Range) -> Bool {
        lowerBound <= otherRange.lowerBound && upperBound >= otherRange.upperBound
    }
}


public class UnitCell {

    var origin:Vector
    var dimensions:[Vector]
    var buffer:[Double]
    var units:[Vector]
    var sizes:[Double]
    var membraneaxis:AXES
    var imagecoords:[Vector]?
    var imageindices:[Int:[Int]]?
    var inverseindices:[Int:Int]?


    public init(_ origin:Vector, _ dimensions:[Vector], _ buffer:[Double], _ membraneaxis:AXES )  {
        self.origin = origin
        self.dimensions = dimensions
        // note that unit vector is an optional
        self.units = dimensions .map { $0.unit()! }
        self.sizes = dimensions .map { $0.length() }
        // translate buffer distances into cell coordinates
        self.buffer = zip( buffer, self.sizes ) .map { $0 / $1 }
        self.membraneaxis = membraneaxis

        self.imagecoords = nil
        self.imageindices = nil
        self.inverseindices = nil

    }

    public func cellcoords( _ coord:Vector ) -> [Double] {
        
        let disp = coord.sub(self.origin)
        let cellc = (0..<3) .map { disp.dot(self.units[$0]) / self.sizes[$0] }

        return cellc
    }

    public func packcell( _ atompos:[Vector] ) -> [Vector] {

    var transPos = [Vector]()

    // this assumes at most one image of each atom is present

    for coord in atompos {

            let cellcoords = self.cellcoords(coord)
            var disp = coord.sub(self.origin)

            for ax in 0..<3 {
                
                if cellcoords[ax] < 0.0 {
                    disp = disp.add(self.dimensions[ax])
                }

                if cellcoords[ax] > 1.0 {
                    disp = disp.sub(self.dimensions[ax])
                }
            }

            transPos.append(self.origin.add(disp))

        }

        return transPos

    }

    // this adds reflections of atoms that are within 'buffer' cell coordinate of the cell boundaries
    // in the 'non-cell membrane' directions
    // 
    // a coordinate can have one image if it is within one buffer boundary, or three if it is close to two boundaries

    public func buildimages( _ coords:[Vector] ) -> ([Vector], [Int:[Int]], [Int:Int]) {

        var bcoords = [Vector]()
        var bmap = [Int:[Int]]()
        var binvmap = [Int:Int]()

        

        for (aidx,coord) in coords.enumerated() {
            let cellc = self.cellcoords( coord )
            var disps = [Vector]()
    

            for ax in 0..<3 {
                if ax == self.membraneaxis.rawValue {
                    continue
                }

                if cellc[ax] < self.buffer[ax] {
                    disps.append(self.dimensions[ax].scale(1.0)) 

                }

                if 1.0 - cellc[ax] < self.buffer[ax] {
                    disps.append(self.dimensions[ax].scale(-1.0))
                }
            }

            if disps.count > 0  {
                if bmap[aidx] == nil {
                    bmap[aidx] = [Int]()
                }
                // apply displacements single and in combination
                var apply = [Vector]()
                if disps.count == 2 {
                    apply.append(disps[0])
                    apply.append(disps[1])
                    apply.append(disps[0].add(disps[1]))
                }
                else {
                    apply.append(disps[0])
                }

                for disp in apply {
                    bmap[aidx]!.append(bcoords.count)
                    bcoords.append(coord.add(disp))
                }
                
            }
        }

        // for convenience make inverse map (image index to original atom index)

        for aidx in bmap.keys {
            for imgidx in bmap[aidx]! {
                binvmap[imgidx] = aidx
            }
        }
        return (bcoords, bmap, binvmap)

    }

    public func setimages( _ imagecoords:[Vector], _ imageindices:[Int:[Int]], _ inverseindices:[Int:Int]) {
        self.imagecoords = imagecoords
        self.imageindices = imageindices
        self.inverseindices = inverseindices
    }

}



public func membraneCoordinates(_ incoordinates:[Vector], _ inradii:[Double], _ proberad:Double, _ unitcell:UnitCell ) 
    -> ([Vector], ([Vector], [Int:[Int]], [Int:Int]), [Double])  {

        var packedcoords = unitcell.packcell( incoordinates )

        let imgdata = unitcell.buildimages(packedcoords) 

        unitcell.setimages( imgdata.0, imgdata.1, imgdata.2 )

        var imgwithidx = [(Int,Int)]()

        for (aidx,imgindices) in unitcell.imageindices! {
            for imgidx in imgindices {
                imgwithidx.append( (imgidx,aidx ))
            }
        }

        let imginorder = imgwithidx .sorted { $0.0 < $1.0 }

        let imgradii = imginorder .map { inradii[$0.1] }

        return (packedcoords, imgdata, imgradii)

}

public func processMembraneProbes( _ probes:[Probe], _ proberad:Double, _ unitcell:UnitCell ) 
            -> ([Probe], [Probe], [Probe]) {

    // in first step, remove probes that lie outside the unit cell 

    var keepprobes = [Probe]()

    for probe in probes {
        let cellc = unitcell.cellcoords(probe.center)
        var outside = false
        for ax in 0..<3 {
            if ax == unitcell.membraneaxis.rawValue {
                continue
            }
            if cellc[ax] > 1.0 || cellc[ax] < 0.0 {
                outside = true
            }
        }
        if outside {
            continue
        }

        keepprobes.append(probe)

    }

    // second step is to generate image probes using the unit cell buffer dimensions

    let probecoords = keepprobes .map { $0.center }

    let probeimgdata = unitcell.buildimages( probecoords )

    let probeinvimg = probeimgdata.2 

    var bufferPROBES = [Probe]() 

    for imgidx in probeinvimg.keys {
        let probidx = probeinvimg[imgidx]!
        let probe = keepprobes[probidx]
        let bufprobe = Probe(center:probeimgdata.0[imgidx], radius:probe.proberadius, atoms:probe.atoms, 
            singleton:probe.singleton, clockwise:probe.clockwisecontour)
        bufferPROBES.append(bufprobe)


    }

    let totPROBES = keepprobes + bufferPROBES

    return (totPROBES, keepprobes, bufferPROBES)


}

public func faceNormal( _ faceidx:Int, _ VERTICES:[Vector], _ FACES:[[Int]], threshold:Double=1.0e-8 ) -> Vector? {

    let verts = (0..<3) .map { FACES[faceidx][$0] }
    let coords = (0..<3) .map { VERTICES[verts[$0]] }

    let edge0 = coords[1].sub(coords[0])
    let edge1 = coords[2].sub(coords[0])

    let orient = edge0.cross(edge1)

    if orient.length() < threshold {
        return nil ;
    }

    return orient.unit()

}

public func projectIntoElementPlane( _ point:Vector,  _ faceidx:Int, _ VERTICES:[Vector], _ FACES:[[Int]],
 _ unitcell:UnitCell) -> Vector? {

    let norm = faceNormal( faceidx, VERTICES, FACES )

    let axis = unitcell.membraneaxis.rawValue
    let RIGHT = axisRIGHT[axis]
    let UP =    axisUP[axis]

    if norm == nil || norm!.coords[axis] == 0.0 {
        return nil
    }

    let v0 = VERTICES[FACES[faceidx][0]]

    let pu = point.coords[RIGHT]
    let pv = point.coords[UP]
    let v0u = v0.coords[RIGHT]
    let v0v = v0.coords[UP]

    let paxis = v0.coords[axis] - ( (pu - v0u)*norm!.coords[RIGHT] + (pv - v0v)*norm!.coords[UP] )/norm!.coords[axis]

    var proj = Vector([0.0,0.0,0.0])

    proj.coords[RIGHT] = pu 
    proj.coords[UP] = pv 
    proj.coords[axis] = paxis

    return proj

}

// assume that unit cell UP limits are high enough that no 'horizontal' axis meets any elements

public func cornerInsideElement( _ corner:Vector, _ faceidx:Int, _ VERTICES:[Vector], _ FACES:[[Int]],
            _ unitcell:UnitCell ) -> Bool {

    let point = projectIntoElementPlane( corner, faceidx, VERTICES, FACES, unitcell )

    if point == nil {
        return false 
    }

    let orient = faceNormal( faceidx, VERTICES, FACES )

    if orient == nil {
        return false
    }

    var violate = false

    let vertices = (0..<3) .map { FACES[faceidx][$0] }
 
    for edgespec in [(0,1),(1,2),(2,0)] {
        let v0 = vertices[edgespec.0]
        let v1 = vertices[edgespec.1]
        let t = point!.sub(VERTICES[v0])
        let e = VERTICES[v1].sub(VERTICES[v0])
        if e.cross(t).dot(orient!) <= 0 {
            violate = true 
            break
        }
    }

    return !violate
}

public func coord3DToUVW( _ point:Vector, _ unitcell:UnitCell ) -> Vector {
    let delta = point.sub(unitcell.origin)
    let u = delta.coords[0]/unitcell.sizes[0]
    let v = delta.coords[1]/unitcell.sizes[1]
    let w = delta.coords[2]/unitcell.sizes[2]

    return Vector([u,v,w])
}

public func uvwToCoord3D( _ point:Vector, _ unitcell:UnitCell ) -> Vector {
    let coord = unitcell.origin.add(unitcell.dimensions[0].scale(point.coords[0]))
        .add(unitcell.dimensions[1].scale(point.coords[1]))
        .add(unitcell.dimensions[2].scale(point.coords[2]))
    
    return coord

}



public func edgeIntersections( _ edge:(Int,Int), _ uvwVERTICES:[Vector] ) -> ([Double],[Vector]) {

    let v0 = edge.0 
    let v1 = edge.1

    var data:([Double],[Vector]) = ([],[])

    var Ts = [Double]()
    var Axs = [Int]()

    // just check all possibilities - 
    // remember that two intersections for same edge may involve different axes or 0/1 crossing, so 
    // cannot prematurely break

    // Also, possible for a value of t to satisfy one unit cell boundary and exceed another

    for ax in 0..<3 {

        let x0 = uvwVERTICES[v0].coords[ax]
        let x1 = uvwVERTICES[v1].coords[ax]

        if x0 == x1 {
            continue
        }

        // 0 cross 

        var t = -x0/(x1 - x0)

        if 0.0 < t && t < 1.0 {
            Ts.append(t)
            Axs.append(ax)
        }

        // 1 cross 

        t = (1 - x0)/(x1 - x0)

        if 0.0 < t && t < 1.0 {
            Ts.append(t)
            Axs.append(ax)
        }
    }

    for (ax,t) in zip(Axs,Ts) {

        var fail = false

        for ax2 in 0..<3 {
            if ax2 == ax { continue }
            let y = (1.0 - t)*uvwVERTICES[v0].coords[ax2] + t*uvwVERTICES[v1].coords[ax2]
            if y < 0.0 || y > 1.0 {
                fail = true 
                break
            }
        }

        if !fail {
                data.0.append(t)
                data.1.append(uvwVERTICES[v0].scale(1.0-t).add(uvwVERTICES[v1].scale(t)))
        }
    }


    return data 
}

public func checkForTransverse( _ faceidx:Int, _ uvwVERTICES:[Vector], _ FACES:[[Int]],
            _ unitcell:UnitCell ) -> ([Int],[Range<Double>]) {
    
    let axis = unitcell.membraneaxis.rawValue 
    let RIGHT = axisRIGHT[axis]
    let UP = axisUP[axis]

    var data:([Int],[Range<Double>]) = ([],[])

    let verts = (0..<3) .map { FACES[faceidx][$0] }

    for (eidx,edge) in [(0,1),(1,2),(2,0)].enumerated() {
        let v0 = verts[edge.0]
        let v1 = verts[edge.1]

        let inters = edgeIntersections( (v0,v1), uvwVERTICES )

        if inters.0.count == 2 {
            let ts = inters.0.sorted { $0 < $1 }
            if ts[0] != ts[1] {
                data.0.append(eidx)
                data.1.append(ts[0]..<ts[1])
            }
        }
    }


    return data
}


public func processMembraneTri( _ VERTICES_in:[Vector], _ NORMALS_in:[Vector], _ FACES_in:[[Int]], _ unitcell:UnitCell ) 
        -> ([Vector], [Vector], [[Int]]) {

    // make mutating 

    var VERTICES = Array(VERTICES_in)
    var NORMALS =  Array(NORMALS_in)
    var FACES =    Array(FACES_in)

    // get corners of unit cell, we only want the implied vertical axes 
    // order : ll, lr, ul, ur

    let axis = unitcell.membraneaxis.rawValue
    let RIGHT = axisRIGHT[axis]
    let UP = axisUP[axis]

    // Note that I am replicating here - I have a cellcoord function in UnitCell, but it returns [Double], 
    // I prefer a Vector here. Will refactor later

    var cornerLL = Vector([0.0,0.0,0.0])
    var cornerLR = Vector([0.0,0.0,0.0])
    var cornerUL = Vector([0.0,0.0,0.0])
    var cornerUR = Vector([0.0,0.0,0.0])

    cornerLL.coords[UP] = 0.0 
    cornerLR.coords[UP] = 0.0
    cornerUL.coords[UP] = 1.0
    cornerUR.coords[UP] = 1.0
    cornerLL.coords[RIGHT] = 0.0
    cornerUL.coords[RIGHT] = 0.0
    cornerLR.coords[RIGHT] = 1.0
    cornerUR.coords[RIGHT] = 1.0

    let corners3D = [ cornerLL, cornerLR, cornerUL, cornerUR ] .map { uvwToCoord3D( $0, unitcell ) }
    let cornerUVW = [ cornerLL, cornerLR, cornerUL, cornerUR ]


    // Here is an approach - start working in 'u,v' coordinates, which are just cell coordinates in the X-Y plane
    // 
    // Edges map to new vertices they produce owing to intersections ; vertices are also produced in the rare 
    // instance that a corner axis of the unit cell passes through an element
    // 
    // Terminology :
    //      inelems === totally inside unit cell 
    //      outelems === element area does not intersect unit cell 
    //      crosselems === elements with vertices both inside and outside unit cell 
    //      cornerelems === elements that enclose a corner axis of the unit cell 
    //      transelems === elements with all vertices outside unit cell which nontheless have intersection 
    //                          with the unit cell 
    // 
    //      To detect 'transverse' crossing :
    //          first try :
    //              for each edge,
    //              u = (1 - t)*u0 + t*u1 = u0 + t*(u1 - u0)
    //              0 < u ==> 
    //                      t > -u0/(u1 - u0) | u1 > u0  -> rangeU_0 = ( -u0/(u1 - u0), inf )
    //                      t < -u0/(u1 - u0) | u1 < u0  -> rangeU_0 = ( -inf, -u0/(u1 - u0))
    //
    //              u < 1 ==>
    //                      t < (1 - u0)/(u1 - u0) | u1 > u0 -> rangeU_1 = ( -inf, (1 - u0)/(u1 - u0))
    //                      t > (1 - u0)/(u1 - u0) | u1 < u0 -> rangeU_1 = ( (1 - u0)/(u1 - u0), inf )
    //              So, for an edge we get four ranges, treating separately u and v ; finaly, the intersection of those
    //              (if it exists) must overlap (0,1)
    //
    //          second try (simpler) :
    //              for each edge :
    //                  find all intersections
    //                  identify edges with two distinct intersections, expect two such edges
    //      If the elements does not contain a corner, logically both edges must be transverse, or both not transverse

    var inelems = [Int]()
    var outelems = [Int]()
    var crosselems = [ (Int, [Int]) ]()   // retain [fidx, inside vertices]
    var cornerelems = [(Int, Int, [Int])]()  // retain [fidx, corner index, inside vertices]
    var transelems:[(Int,([Int],[Range<Double>]))] = []   // retain [fidx, <intersection data>]
                                                          // intersection data is pairs, edge position (0,1,2) paired with range of t for intersections
                                                          // parameter t corresponds to original edge orientation 
    
    
    let vertexCellCoords = VERTICES .map { coord3DToUVW( $0, unitcell ) }

    // Note that an element can overlap the interior of the X-Y cell even if all vertices lie outside, need 
    // additional test

    for (fidx,face) in FACES.enumerated() {
        var inside = [Int]()
        var outside = [Int]()
        var minU =  Double.infinity
        var maxU = -Double.infinity
        var minV =  Double.infinity 
        var maxV = -Double.infinity

        for vidx in face  {

            let vertex = vertexCellCoords[vidx]

            var out = false

            minU = [minU, vertex.coords[RIGHT]].min()!
            maxU = [maxU, vertex.coords[RIGHT]].max()!
            minV = [minV, vertex.coords[UP]].min()!
            maxV = [maxV, vertex.coords[UP]].max()!
         

            for ax in 0..<3 {
                if ax == unitcell.membraneaxis.rawValue { continue }
                if vertex.coords[ax] < 0.0 || vertex.coords[ax] > 1.0 {
                    out = true
                    break
                }
            }
            
            if out {
                outside.append(vidx)
            }
            else {
                inside.append(vidx)
            }

        } 


        if outside.count == 0 {
            inelems.append(fidx)
        }
        else {

            // check if possibility of corner containment 

            var cornerinside = false

            for (cidx,corner) in cornerUVW.enumerated() {
                if minU < corner.coords[RIGHT] && minV < corner.coords[UP] && maxU > corner.coords[RIGHT] 
                    && maxV > corner.coords[UP] {
                    // check for interior 
                    cornerinside = cornerInsideElement( corners3D[cidx], fidx, VERTICES, FACES, unitcell )
                    if cornerinside {
                        cornerelems.append((fidx, cidx, inside))
                        break
                    }
                }
            }

            if !cornerinside {
                if inside.count == 0 {
                    let tdata = checkForTransverse( fidx, vertexCellCoords, FACES,
                             unitcell )
                    if tdata.0.count == 2 {
                        transelems.append((fidx,tdata))
                    }
                    else if tdata.0.count == 0 {
                        outelems.append(fidx)
                    }
                    else {
                        print("\nwarning, unexpected : element \(fidx) has \(tdata.0.count) transverse intersections")
                    }
                    
                }
                else {
                    crosselems.append((fidx, inside))
                }
            }
        }

    }

    // map edges to potential new vertices generated by intersections 
    // 
    // sort edge vertices as we need a unique identifier

    // edge maps to t, 3d intersection point, vertex assignment (-1 initially), and list of face indices contributing to intersection

    var edgeToIntersections3D = [Edge:[(Double,Vector,Int,[Int])]]()

    // handle each element type with intersections separately ; for each one, map edges to 3D intersections, and 
    // value of parameter t for edge vertices in **sorted order**

    // combine all element indices where we need to find edge intersections

    var interelems = [Int]()

    crosselems  .map { interelems.append( $0.0 ) }
    cornerelems .map { interelems.append( $0.0 ) }
    transelems  .map { interelems.append( $0.0 ) }

    // note that edge intersections are only found with the unit cell boundaries;
    // so, it is enough to check each edge only once  

    var edgeChecked = Set<Edge>()

    for fidx in interelems {

        for edgespec in [(0,1),(1,2),(2,0)] {
            
            let v0 = FACES[fidx][edgespec.0]
            let v1 = FACES[fidx][edgespec.1]
            let sortverts = [v0,v1].sorted { $0 < $1 }
            let theEdge = Edge(v1:sortverts[0], v2:sortverts[1])

            if edgeChecked.contains(theEdge) {

                // only need to add fidx to contributing elements 

                if edgeToIntersections3D[theEdge] != nil {
                    for k in 0..<edgeToIntersections3D[theEdge]!.count {
                        edgeToIntersections3D[theEdge]![k].3.append(fidx)
                    }
                }

                continue 
            }

            edgeChecked.insert(theEdge)

            var inters = edgeIntersections( (v0,v1), vertexCellCoords )

            if inters.0.count > 0 {
                
                let reverse = sortverts[0] != v0
                var interpoints = [(Double,Vector,Int,[Int])]()
                // it is enough to recompute parameter t here to correspond to edge reversal;
                // we have 3D points computed from t, and eventually all intersections are sorted by t along each edge
                if reverse {
                    for (tidx,t) in inters.0.enumerated() {
                        inters.0[tidx] = 1.0 - t
                    }
                }
                for t in inters.0 {
                    let p = VERTICES[sortverts[0]].scale(1.0-t).add(VERTICES[sortverts[1]].scale(t))
                    interpoints.append((t,p,-1,[fidx]))
                }

                edgeToIntersections3D[theEdge] = interpoints
                
            }
        }

    }

    // all intersections generated, sort them by t along each edge , add corresponding vertices

    // need to estimate normal for added vertices


    for edge in edgeToIntersections3D.keys {

        edgeToIntersections3D[edge]! = edgeToIntersections3D[edge]!.sorted { $0.0 < $1.0 }

        for k in 0..<edgeToIntersections3D[edge]!.count {

            // changing a COPY
            //var inter = edgeToIntersections3D[edge]![k]
            let vidx = VERTICES.count 
            edgeToIntersections3D[edge]![k].2 = vidx 
            VERTICES.append(edgeToIntersections3D[edge]![k].1)

            // estimate normal 

            var normalsum = Vector([0.0,0.0,0.0])

            var count = 0 

            for fidx in edgeToIntersections3D[edge]![k].3 {

                let fn = faceNormal( fidx, VERTICES, FACES )

                if fn == nil {
                    continue
                }

                count += 1
                normalsum = normalsum.add(fn!)

            }

            if count == 0 {
                print("\nwarning : could not estimate normal for new vertex")
                continue
            }

            NORMALS.append( normalsum.unit()! )
        }
    }

    // have all intersections, need special handling for element types 

    var removedElements = [Int]()


    for data in cornerelems {

        let fidx = data.0 
        let cidx = data.1

        let center = projectIntoElementPlane( corners3D[cidx], fidx, VERTICES, FACES, unitcell )
        
        let vidx = VERTICES.count 

        VERTICES.append(center!)

        // take normal as the face normal

        let fn = faceNormal( fidx, VERTICES, FACES )

        NORMALS.append(fn)

        // expect total of two intersections ;

        var intersections:[(Double,Vector,Int,[Int])] = []

        for edgespec in [(0,1),(1,2),(2,0)] {
            let verts = [ FACES[fidx][edgespec.0], FACES[fidx][edgespec.1]]
            let sortedverts = verts.sorted { $0 < $1 }
            let reverse = verts[0] != sortedverts[0]
            let edge = Edge(v1:sortedverts[0], v2:sortedverts[1])
            if edgeToIntersections3D[edge] != nil {
                if reverse {
                    intersections += edgeToIntersections3D[edge]!.reversed()
                }
                else {
                    intersections += edgeToIntersections3D[edge]!
                }
            }
        } 

        if intersections.count != 2 {
            print("\nerror in corner element processing, have \(intersections.count) edge intersections")
            continue
        }

        FACES.append([vidx,intersections[0].2, intersections[1].2] )

        removedElements.append(fidx)

    }


    for data in transelems {

        let fidx = data.0 

        var intersections:[(Double,Vector,Int,[Int])] = []

        for edgespec in [(0,1),(1,2),(2,0)] {
            let verts = [ FACES[fidx][edgespec.0], FACES[fidx][edgespec.1]]
            let sortedverts = verts.sorted { $0 < $1 }
            let reverse = verts[0] != sortedverts[0]
            let edge = Edge(v1:sortedverts[0], v2:sortedverts[1])
            if edgeToIntersections3D[edge] != nil {
                if reverse {
                    intersections += edgeToIntersections3D[edge]!.reversed()
                }
                else {
                    intersections += edgeToIntersections3D[edge]!
                }
            }
        } 

        if intersections.count != 4 {
            print("\nerror in transverse element processing, have \(intersections.count) edge intersections")
            continue
        }

        // connect 0-1-2, 2-3-0

        FACES.append([intersections[0].2, intersections[1].2, intersections[2].2] )
        FACES.append([intersections[2].2, intersections[3].2, intersections[0].2] )

        removedElements.append(fidx)

    }

    for data in crosselems {

        // this is the most complex case, I believe

        // can have one or two vertices inside

        // order edges based on number of intersections

        var intersectionsForEdges:[[(Double,Vector,Int,[Int])]?] = [nil, nil, nil]

        let fidx = data.0 

        var insideInOrder = [Int]()

        if data.1.count == 1 {
            insideInOrder.append(data.1[0])
        }

        for (eidx,edgespec) in ([(0,1),(1,2),(2,0)]).enumerated() {
            let verts = [ FACES[fidx][edgespec.0], FACES[fidx][edgespec.1]]
            if data.1.count == 2 {
                if data.1.contains(verts[0]) && data.1.contains(verts[1]) {
                    insideInOrder += verts
                }
            } 
            let sortedverts = verts.sorted { $0 < $1 }
            let reverse = verts[0] != sortedverts[0]
            let edge = Edge(v1:sortedverts[0], v2:sortedverts[1])
            if edgeToIntersections3D[edge] != nil {
                if reverse {
                    intersectionsForEdges[eidx] = Array(edgeToIntersections3D[edge]!.reversed())
                }
                else {
                    intersectionsForEdges[eidx] = edgeToIntersections3D[edge]!
                }
            }
        } 

        // Note there is an unusual case near cell corner where one vertex is inside and two are outside

        let intersections = intersectionsForEdges .filter { $0 != nil }

        if intersections.count != 2 &&  intersections.count != 4 {
            print("\nerror in crossing element processing elem # \(fidx), have \(intersections.count) edge intersections")
            continue
        }


        //   if any edge has two intersections, rotate so it is in position 0
        // ELSE 
        //   if two interior points, rotate so that 'no intersection' edge is in position 2
        //      ELSE
        //   move 'no intersection' edge to position 1

        var twointersection:Int?
        var nointersection:Int? 

        for j in 0..<3 {
            if intersectionsForEdges[j] != nil {
                if intersectionsForEdges[j]!.count == 2 {
                    twointersection = j
                }
            }
            else {
                nointersection = j
            }
        }

        var rotate = 0

        if twointersection != nil {
            rotate = 3 - twointersection!
        }
        else {
            if insideInOrder.count == 1 {
                // nointersection to position 1
                if nointersection! == 0 {
                    rotate = 1
                }
                else if nointersection! == 2 {
                    rotate = 2
                }
            }
            else {
                // nointersection to position 2\

                rotate = 2 - nointersection!

            }
            
        }

        var edgespec = [(0,1),(1,2),(2,0)]

        if rotate == 1 {
            edgespec = [ edgespec[2], edgespec[0], edgespec[1]]
            intersectionsForEdges = [intersectionsForEdges[2], intersectionsForEdges[0], intersectionsForEdges[1] ]
        }
        else if rotate == 2 {
            edgespec = [ edgespec[1], edgespec[2], edgespec[0]]
            intersectionsForEdges = [intersectionsForEdges[1], intersectionsForEdges[2], intersectionsForEdges[0] ]
        }


        
        if data.1.count == 1 {
            if twointersection == nil {
                FACES.append( [insideInOrder[0], intersectionsForEdges[0]![0].2, intersectionsForEdges[2]![0].2] )
            }
            else {
                FACES.append( [ intersectionsForEdges[0]![0].2, intersectionsForEdges[0]![1].2, insideInOrder[0] ] )
                FACES.append( [ intersectionsForEdges[0]![1].2, intersectionsForEdges[1]![0].2, insideInOrder[0] ] )
                FACES.append( [ intersectionsForEdges[2]![0].2, intersectionsForEdges[0]![1].2, insideInOrder[0] ] )
            }
        }
        else {
            
                FACES.append( [insideInOrder[0], intersectionsForEdges[0]![0].2, intersectionsForEdges[1]![0].2 ] )
                FACES.append( [insideInOrder[1], intersectionsForEdges[0]![0].2, insideInOrder[0] ] )
            
        }

        removedElements.append(fidx)

    }



    /*
    for data in crosselems {
        let fidx = data.0
        let insideverts = data.1
        let verts = (0..<3) .map { FACES[fidx][$0] }

        for edgespec in [(0,1),(1,2),(2,0)] {
            let minv = [ verts.edgespec.0,verts.edgespec.1 ].min()!
            let maxv = [ verts.edgespec.0,verts.edgespec.1 ].max()!
            if !insideverts.contains(minv) && !insideverts.contains(maxv) { continue }
            let reverse = verts[edgespec.1] != maxv
            let edge = Edge(v1:minv, v2:maxv)

            let uc0 = uvwVERTICES[minv]
            let uc1 = uvwVERTICES[minv]

            var t = -1.0

            for ax in [RIGHT, UP]  {
                let x0 = uc0.coords[ax]
                let x1 = uc1.coords[ax]
                if x0*x1 < 0.0 {
                    // have left or bottom crossing
                    t = (-x0)/(x1 - x0)
                }
                else if (1.0 - x0)*(1.0 - x1) < 0.0 {
                    // have right or top crossing
                    t = (1.0 - x0)/(x1 - x0)
                }
                else {
                    print("\ncannot determine crossing type for edge connecting vertices \(minv) - \(maxv)")
                    print("\tuc coords \(minv) = \(uc0.coords[0]),\(uc0.coords[1]),\(uc0.coords[2])")
                    print("\tuc coords \(maxv) = \(uc1.coords[0]),\(uc1.coords[1]),\(uc1.coords[2])")
                }
            }
            if t < 0 { continue }

            let intersection = VERTICES[minv].scale(1.0 - t).add(VERTICES[maxv].scale(t))

            if edgeToIntersections[edge] == nil {
                edgeToIntersections[edge] = [ (intersection,reverse,t) ]
            }
            else {
                edgeToIntersections[edge]!.append( (intersection,reverse,t) )
            }
        }
    }

    */


    // find intersections for corner

    // not sure that initial approach of simply deleting elements at one side is adequate - what happens upon decimation??
    // 
    // instead add edges to make boundary elements conform to unit cell boundaries; means adding vertices, edges and faces

    // for oriented element [0, 1, 2, 0] :
    // if two vertices separated by u.c. boundary, insert a vertex between 
    // e.g. if 2 and 0 on opp sides, have [0, 1, 2, 20, 0 ]
    // have faces [0, 1, 2] + [0, 2, 20]
    //
    // if 0, 2 opp 1 : [0, 01, 1, 2, 20, 0], elements [0, 01, 20], [01, 1, 20], [1, 2, 20]
    //
    // special condition - does face contain unit cell 'corner' ? This gets really complicated; maybe see first if 
    // smoothing and decimation looks OK with the ragged edge.

    /*
    var newfaces = [[Int]]()

    for fidx in crosselems {

        let insideverts =  FACES[fidx] .filter { vertexcellcoords[$0] <= 1.0 || vertexcellcoords[$0] >= 0.0 }
        let outsideverts = FACES[fidx] .filter { vertexcellcoords[$0] > 1.0 || vertexcellcoords[$0] < 0.0 }

        var extendedface = [Int]()

        for i in 0..<3 {
            let j = i < 2 ? i + 1 : 0

            let vi = FACES[fidx][i]
            let vj = FACES[fidx][j]

            if (insideverts.contains(vi) && outsideverts.contains(vj)) || (insideverts.contains(vj) && outsideverts.contains(vi)) {
                
            }
        }



    } 

    */


    // see if I can identify symmetry-related elements 

    // for convenience, get centroids in uc coords for the elements
    /*

    let centroids = FACES .map { vertexCellCoords[$0[0]].add(vertexCellCoords[$0[1]]).add(vertexCellCoords[$0[2]]).scale( 1.0/3.0) }

    var OUT = "\n\n*** CORNER ELEMENTS (unit cell coordinates)\n"

    cornerelems = cornerelems.sorted { (centroids[$0.0].coords[2],  centroids[$0.0].coords[0], centroids[$0.0].coords[1]) < 
            (centroids[$1.0].coords[2],  centroids[$1.0].coords[0], centroids[$1.0].coords[1])     }

    for data in cornerelems {
        let verts = (0..<3) .map { FACES[data.0][$0] }
        var outstr1 = "\n\t verts: "
        var outstr2 = "\t u,v,w : "
        for v in verts {
            outstr1 += " \(v)"
            outstr2 += " \(vertexCellCoords[v].coords[0]) \(vertexCellCoords[v].coords[1]) \(vertexCellCoords[v].coords[2])"
        }
        outstr1 += "\n"
        outstr2 += "\n"
        OUT += outstr1
        OUT += outstr2
    }

     OUT += "\n\n*** CROSS ELEMENTS (unit cell coordinates)\n"

    crosselems = crosselems.sorted { (centroids[$0.0].coords[2],  centroids[$0.0].coords[0], centroids[$0.0].coords[1]) < 
            (centroids[$1.0].coords[2],  centroids[$1.0].coords[0], centroids[$1.0].coords[1])     }

    for data in crosselems {
        let verts = (0..<3) .map { FACES[data.0][$0] }
        var outstr1 = "\n\t verts: "
        var outstr2 = "\t u,v,w : "
        for v in verts {
            outstr1 += " \(v)"
            outstr2 += " \(vertexCellCoords[v].coords[0]) \(vertexCellCoords[v].coords[1]) \(vertexCellCoords[v].coords[2])"
        }
        outstr1 += "\n"
        outstr2 += "\n"
        OUT += outstr1
        OUT += outstr2
    }

     OUT += "\n\n*** TRANSVERSE ELEMENTS (unit cell coordinates)\n"

    transelems = transelems.sorted { (centroids[$0.0].coords[2],  centroids[$0.0].coords[0], centroids[$0.0].coords[1]) < 
            (centroids[$1.0].coords[2],  centroids[$1.0].coords[0], centroids[$1.0].coords[1])     }

    for data in transelems {
        let verts = (0..<3) .map { FACES[data.0][$0] }
        var outstr1 = "\n\t verts: "
        var outstr2 = "\t u,v,w : "
        for v in verts {
            outstr1 += " \(v)"
            outstr2 += " \(vertexCellCoords[v].coords[0]) \(vertexCellCoords[v].coords[1]) \(vertexCellCoords[v].coords[2])"
        }
        outstr1 += "\n"
        outstr2 += "\n"
        OUT += outstr1
        OUT += outstr2
    }

    let url = URL(fileURLWithPath: "./element_output.txt")

    do {
        try OUT.write(to: url, atomically: true, encoding: String.Encoding.utf8)
    } catch {
        print("error writing file \(url)")
    }


    // delete crossover elements at 'right' or 'top' of unit cell --- 
    // in brief, keep any element with all vertex coordinates less than 1.0

    /*
    var keepelems = Array(inelems) 

    for data in crosselems {
        let fidx = data.0
        var keep = true 
        let vertcoord = FACES[fidx] .map { vertexCellCoords[$0] }
        
        for (aidx,coord) in vertcoord.enumerated() {
            if aidx == unitcell.membraneaxis.rawValue { continue }
            if coord.coords[aidx] > 1.0 {
                keep = false
                break
            }
        }

        if keep {
            keepelems.append(fidx)
        }
    }

    var vertexindices = [Int]()

    for fidx in keepelems {
        vertexindices += FACES[fidx]
    }

    let keepvertexset = Set(vertexindices)

    var keepvertexindices = Array(keepvertexset) 
    keepvertexindices = keepvertexindices .sorted { $0 < $1 }

    // translate old to new vertex indices 

    var translate = Array( repeating:-1 , count:VERTICES.count )

    for (newidx,oldidx) in keepvertexindices.enumerated() {
        translate[oldidx] = newidx
    }

    // vertex data for export 

    let exportVertices = keepvertexindices .map { VERTICES[$0] }
    let exportNormals =  keepvertexindices .map { NORMALS[$0] }

    var exportElems = [[Int]]()

    for fidx in keepelems {
        let v0 = FACES[fidx][0]
        let v1 = FACES[fidx][1]
        let v2 = FACES[fidx][2]
        exportElems.append( [translate[v0], translate[v1], translate[v2]])
    }
    */
    */

    // keep elems are :
    // any not in removedElements
    // any with index >= ORIG_FACE_COUNT
    // any not in outelems

    var keepElems = Array( repeating:true, count:FACES.count )

    for fidx in removedElements {
        keepElems[fidx] = false 
    }

    for fidx in outelems {
        keepElems[fidx] = false
    }

    var keepVertices = Array( repeating:false, count:VERTICES.count )

    for (fidx,keep) in keepElems.enumerated() {
        if !keep {
            continue
        }

        for v in FACES[fidx] {
            keepVertices[v] = true
        }
    }

    var vertexToNewIndex = Array( repeating:-1, count:VERTICES.count )

    var currentIndex = 0 

    for (vidx,keep ) in keepVertices.enumerated() {
        if !keep { 
            continue
        }

        vertexToNewIndex[vidx] = currentIndex 
        currentIndex += 1
    }

    let exportVertices = VERTICES.enumerated() .filter { vertexToNewIndex[$0.offset] >= 0  } .map { $0.element } 
    let exportNormals = NORMALS.enumerated() .filter { vertexToNewIndex[$0.offset] >= 0  } .map { $0.element }
    let exportElems = FACES .filter { vertexToNewIndex[$0[0]] >= 0 && vertexToNewIndex[$0[1]] >= 0 && vertexToNewIndex[$0[2]] >= 0 }
            .map { [vertexToNewIndex[$0[0]], vertexToNewIndex[$0[1]], vertexToNewIndex[$0[2]] ]}


    return (exportVertices, exportNormals, exportElems)


}



public func mergeSurfaceComponents( _ vertA:[Vector], _ vertB:[Vector], _ normA:[Vector], _ normB:[Vector],
    _ faceA:[[Int]], _ faceB:[[Int]]) -> ([Vector],[Vector],[[Int]]) {

        var outverts = [Vector]() 

        var outnorms = [Vector]()

        outverts += vertA 
        outverts += vertB

        outnorms += normA 
        outnorms += normB 

        let vertexmapB = (0..<vertB.count) .map { $0 + vertA.count }

        var outfaces = [[Int]]()

        outfaces += faceA 

        outfaces += faceB .map { [vertexmapB[$0[0]], vertexmapB[$0[1]], vertexmapB[$0[2]]] }

        return ( outverts, outnorms, outfaces ) 
    } 



public func membraneSurfaceComponents( _ SUBVERTICES:[[Vector]], _ SUBNORMALS:[[Vector]], _ SUBFACES:[[[Int]]], 
            _ unitcell:UnitCell)
             ->  ([[Vector]], [[Vector]], [[[Int]]])  {

        // strategy - take random sample of biggest components, for a sample of elements get normal along membrane axis, 
        // average coordinate in direction of axis; expect highest and lowest with opposite normals, inner two with likwise 
        // opposite normals ; merge outer and inner into single components ('0' is outer, '1' is inner)

        var OUTVERTICES = [[Vector]]()
        var OUTNORMALS = [[Vector]]()
        var OUTFACES = [[[Int]]]()

        var maxPositions = [Double]()

        let ax = unitcell.membraneaxis.rawValue
        

        for isurf in 0..<4 {
            

            var axcoors = SUBVERTICES[isurf] .map { $0.coords[ax] }
            
            maxPositions.append( axcoors.max()! )

        }

        

        var compAxisPos = (0..<4) .map { (maxPositions[$0],$0) } 

        compAxisPos = compAxisPos .sorted { $0.0 < $1.0 }

        let bottom0 = compAxisPos[0].1
        let bottom1 = compAxisPos[1].1
        let top1 = compAxisPos[2].1
        let top0 = compAxisPos[3].1

        
        // membrane components are sorted - merge top.0 and bottom.0, top.1 and bottom.1 ; any
        // other components are copied

        let membrane0 = mergeSurfaceComponents( SUBVERTICES[top0], SUBVERTICES[bottom0], 
                SUBNORMALS[top0], SUBNORMALS[bottom0], SUBFACES[top0], SUBFACES[bottom0] )

        let membrane1 = mergeSurfaceComponents( SUBVERTICES[top1], SUBVERTICES[bottom1], 
                SUBNORMALS[top1], SUBNORMALS[bottom1], SUBFACES[top1], SUBFACES[bottom1] )

        OUTVERTICES.append(membrane0.0)
        OUTVERTICES.append(membrane1.0)

        OUTNORMALS.append(membrane0.1)
        OUTNORMALS.append(membrane1.1)

        OUTFACES.append(membrane0.2)
        OUTFACES.append(membrane1.2)

        if SUBVERTICES.count > 4 {
            for isurf in 4..<SUBVERTICES.count {
                OUTVERTICES.append(SUBVERTICES[isurf])
                OUTNORMALS.append(SUBNORMALS[isurf])
                OUTFACES.append(SUBFACES[isurf])
            }
        }

        return (OUTVERTICES,OUTNORMALS,OUTFACES)


    }