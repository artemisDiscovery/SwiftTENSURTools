
import Foundation 
import MathTools

public struct Edge : Hashable {
    var v1:Int
    var v2:Int

    public func hash(into hasher: inout Hasher) {
        hasher.combine(v1)
        hasher.combine(v2)
    }
}

// if we have boundary, cannot really determine topology of structure - 

public enum SurfaceType: Int {
    case probeCenteredClosed    // no boundary, 'outside'
    case reentrantClosed        // no boundary, 'inside'
    case reentrantCavity        // no boundary, 'outside', interior to an 'inside', contains no other surface
    case undeterminedOpen      // has boundary, cannot determined
    
}

// going to insert a fast algorithm for components here, just need faces and total number of vertices (incase some do not form faces)

// return vertices by component, faces by component

// try to approaches, my original recasting of DFS, and the CFS algorithm

public func surfaceComponentsDFS( faces:[[Int]], numvertices:Int ) -> [[Int]] {

    var adjacency = [Set<Int>]()

    for _ in 0..<numvertices {
        adjacency.append(Set<Int>())
    }

    for f in faces {
        adjacency[f[0]].insert(f[1])
        adjacency[f[1]].insert(f[0])
        adjacency[f[0]].insert(f[2])
        adjacency[f[2]].insert(f[0])
        adjacency[f[1]].insert(f[2])
        adjacency[f[2]].insert(f[1])
    }

    var STACK = [[Int]]()

    var visited = Array(repeating:false, count:numvertices)
    var component = Array(repeating:-1, count:numvertices)

    var currentComponent = -1
    var unassigned:Int?

    while true {
        //find first unassigned vertex

        unassigned = nil

        for iv in 0..<numvertices {
            if component[iv] < 0 {
                unassigned = iv
                break
            }
        }
        if unassigned == nil {
            break
        }

        currentComponent += 1

        STACK.append([unassigned!,currentComponent])

        while STACK.count > 0 {
            let data = STACK.popLast()!
            if !visited[data[0]] {
                visited[data[0]] = true 
                component[data[0]] = data[1]
                for iv in adjacency[data[0]] {
                    STACK.append([iv,currentComponent])
                }
            }
        }
    }

    var components = [[Int]]()

    for c in 0..<(currentComponent+1) {
        components.append([Int]())
    }

    for (iv,c) in component.enumerated() {
        components[c].append(iv)
    }

    components = components .sorted { $0.count > $1.count }

    return components

}


public func surfaceComponentsCFS( faces:[[Int]], numvertices:Int ) -> [[Int]] {

    // initialize from edges 

    var adjacency = Array( repeating:[Int](), count:numvertices )

    for face in faces {
        adjacency[face[0]].append(face[1])
        adjacency[face[1]].append(face[0])
        adjacency[face[1]].append(face[2])
        adjacency[face[2]].append(face[1])
        adjacency[face[2]].append(face[0])
        adjacency[face[0]].append(face[2])
    }

    for k in 0..<numvertices {
        adjacency[k] = adjacency[k] .sorted { $0 < $1 }
    }

    var keyValue = [(Int,Int)]()

    var iter = 0

    while iter < 20 {

        var newPairs = 0 

        for k in 0..<numvertices {
            if adjacency[k].count == 0 {
                continue
            }

            let vmin = adjacency[k][0]

            if vmin > k {
                continue
            }

            keyValue.append((k,vmin))

            for v in adjacency[k] {
                if v != vmin {
                    keyValue.append((v,vmin))
                    newPairs += 1
                }
            }


        }

        if newPairs == 0 {
            break
        }

        adjacency = Array( repeating:[Int](), count:numvertices )

        for tup in keyValue {
            adjacency[tup.0].append(tup.1)
            adjacency[tup.1].append(tup.0)
        }

        for k in 0..<numvertices {
            adjacency[k] = Array( Set(adjacency[k]) ) .sorted { $0 < $1 }
        }
        
    }

    // assemble components 

    var components = [Int:[Int]]() 

    for pair in keyValue {
        if components[pair.0] == nil {
            components[pair.0] = [ pair.0 ]
        }
        components[pair.0]!.append(pair.1) 
    }

    let roots = components.keys .sorted { components[$0]!.count > components[$1]!.count }

    var COMPONENTS = [[Int]]() 

    for root in roots {
        COMPONENTS.append(components[root]!)
    }

    return COMPONENTS

}


// for identifying 'caps' that sew to main membrane reentrant surfaces 

public func edgeComponentsCFS( boundaryedges:[Edge] ) -> [[Int]] {

    print("\nDEPRECATED : function edgeComponentsCFS : not tested")
    // only have boundary edges - translate vertex indices 

    var indices = [Int]() 

    for edge in boundaryedges {
        indices += [ edge.v1, edge.v2 ]
    }

    indices = Array(Set(indices)) .sorted { $0 < $1 }

    var oldToNewIndex = Array( repeating:-1, count:indices.max()! )
    var newToOldIndex = [Int]()

    for (newidx,oldidx) in indices.enumerated() {
        oldToNewIndex[oldidx] = newidx
        newToOldIndex.append(oldidx)
    }
    
    let numvertices = indices.count 

    let edges = boundaryedges .map { Edge( v1:oldToNewIndex[$0.v1], v2:oldToNewIndex[$0.v2]) }

    // initialize from edges 

    var adjacency = Array( repeating:[Int](), count:numvertices )

    for edge in edges {

        adjacency[edge.v1].append(edge.v2)
        adjacency[edge.v2].append(edge.v1)
    }

    for k in 0..<numvertices {
        adjacency[k] = adjacency[k] .sorted { $0 < $1 }
    }

    var keyValue = [(Int,Int)]()

    var iter = 0

    while iter < 20 {

        var newPairs = 0 

        for k in 0..<numvertices {
            if adjacency[k].count == 0 {
                continue
            }

            let vmin = adjacency[k][0]

            if vmin > k {
                continue
            }

            keyValue.append((k,vmin))

            for v in adjacency[k] {
                if v != vmin {
                    keyValue.append((v,vmin))
                    newPairs += 1
                }
            }

        }

        if newPairs == 0 {
            break
        }

        adjacency = Array( repeating:[Int](), count:numvertices )

        for tup in keyValue {
            adjacency[tup.0].append(tup.1)
            adjacency[tup.1].append(tup.0)
        }

        for k in 0..<numvertices {
            adjacency[k] = Array( Set(adjacency[k]) ) .sorted { $0 < $1 }
        }
        
    }

    // assemble components 

    var components = [Int:[Int]]() 

    for pair in keyValue {
        if components[pair.0] == nil {
            components[pair.0] = [ pair.0 ]
        }
        components[pair.0]!.append(pair.1) 
    }

    let roots = components.keys .sorted { components[$0]!.count > components[$1]!.count }

    var COMPONENTS = [[Int]]() 

    // map back to original index

    for root in roots {
        COMPONENTS.append(components[root]! .map { newToOldIndex[$0] })
    }

    return COMPONENTS

}




public func areaForFace( vertices:[Vector], face:[Int] ) -> (normal:Vector,area:Double) {

    let disp01 = vertices[face[1]].sub(vertices[face[0]])
    let disp02 = vertices[face[2]].sub(vertices[face[0]])

    var normal = disp01.cross(disp02)
    let len = normal.length()
    let area = 0.5 * len

    if area < 0.000001 {
        normal = Vector([0.0,0.0,0.0])
    }
    else {
        normal = normal.scale(1.0/len)
    }

    return (normal:normal, area:area )
}

// this finds faces and renumbers vertices 

public func findComponentFaces( component:[Int], FACES:[[Int]], numvertices:Int ) -> [[Int]] {


    var oldToNewIndex = Array( repeating:-1, count:numvertices )

    for (newidx,oldidx) in component.enumerated() {
        oldToNewIndex[oldidx] = newidx
    }

    var componentFaces = [[Int]]()

    for face in FACES {

        let newids = face .map { oldToNewIndex[$0] }

        if newids.contains(-1) {
            if newids .reduce ( 0, + ) != -3 {
                print("warning, in processComponent mapped element vertices = \(newids)")
            }
            continue
        }

        componentFaces.append(newids)

    }

    return componentFaces

}

public func findBoundaryEdges( faces:[[Int]] ) -> (edges:[Edge],vertices:[Int]) {

    
    var edgeToFaces = [Edge:[Int]]()
    
    for (fidx,face) in faces.enumerated() {

        for edgedesc in [(0,1),(1,2),(2,0)] {
            let verts = [ face[edgedesc.0], face[edgedesc.1] ] .sorted { $0 < $1 }
            let edge = Edge(v1:verts[0], v2:verts[1])
            if edgeToFaces[edge] == nil {
                edgeToFaces[edge] = []
            }
            edgeToFaces[edge]!.append(fidx)
        }
    }

    var boundedges = [Edge]() 

    for edge in edgeToFaces.keys {
        if edgeToFaces[edge]!.count < 2 {
            boundedges.append(edge)
        }
    }

    // return boundary vertices and edges 

    var boundvertices = [Int]()

    for edge in boundedges {
        boundvertices.append( edge.v1 )
        boundvertices.append( edge.v2 )
    }

    boundvertices = Array(Set(boundvertices)) .sorted { $0 < $1 }

    return (edges:boundedges, vertices:boundvertices)
}

public func identifySurfaceType( vertices:[Vector], normals:[Vector], boundaryvertices:[Int]=[], axis:Int=2 ) -> SurfaceType {

    // for now just use Z-axis for simplicity, no longer imagining 'caps'

    if boundaryvertices.count > 0 {
        return SurfaceType.undeterminedOpen
    }

    let maxidx = vertices.indices.max( by: { vertices[$0].coords[axis] < vertices[$1].coords[axis] } )!

    let maxCoord = vertices[maxidx].coords[axis]

    var direction = Vector([0.0,0.0,0.0])

    direction.coords[axis] = 1.0 

    let test = normals[maxidx].dot(direction)

    if test > 0 {
        return SurfaceType.reentrantClosed
    }
    else {
        return SurfaceType.probeCenteredClosed
    }


}


public func reverseNormals( normals:[Vector], faces:[[Int]]) -> (normals:[Vector], faces:[[Int]]) {

    // reverse direction of normals, and flip orientation of faces 

    let reversenormals = normals .map { $0.scale(-1.0) }

    let reversefaces = faces .map { [ $0[0], $0[2], $0[1] ] }

    return ( normals:reversenormals, faces:reversefaces )

}

public func boundingBox( points:[Vector] ) -> (min:Vector,max:Vector) 
 {
    var mins = [Double]()
    var maxs = [Double]()
    
    for ax in 0..<3 {
        let coords = points .map { $0.coords[ax] }
        let min = coords.min()!
        let max = coords.max()!
        mins.append(min)
        maxs.append(max)
    }

    return (min:Vector(mins), max:Vector(maxs))
    
}

public func boundingBoxInside( smaller:(min:Vector,max:Vector), larger:(min:Vector,max:Vector) ) -> Bool {

    for ax in 0..<3 {
        if smaller.min.coords[ax] < larger.min.coords[ax] {
            return false
        }

        if smaller.max.coords[ax] > larger.max.coords[ax] {
            return false
        }
    }

    return true
}

public func isInterior( smaller:(vertices:[Vector],normals:[Vector],faces:[[Int]],surfacetype:SurfaceType), 
    larger:(vertices:[Vector],normals:[Vector],faces:[[Int]],surfacetype:SurfaceType)) -> Bool {

        // bounding box of smaller inside larger ??

        let smallbox = boundingBox( points:smaller.vertices )
        let largebox = boundingBox( points:larger.vertices  )

        if !boundingBoxInside( smaller:smallbox, larger:largebox ) {
            return false 
        }

        // do more precise test, just use Z axis 

        let maxidx = smaller.vertices.indices.max( by: { smaller.vertices[$0].coords[2] < smaller.vertices[$1].coords[2] } )!

        // get top point for smaller surface, convert to bounding box

        let coord = (0..<3) .map { smaller.vertices[maxidx].coords[$0] }

        let maxZ = coord[2]

        let pointbox = (min:Vector([coord[0],coord[1],0.0]),max:Vector([coord[0],coord[1],0.0]))

        var candidateFaces = [Int]()

        for (fidx,face) in larger.faces.enumerated() {

            let points = face .map { Vector([larger.vertices[$0].coords[0],  larger.vertices[$0].coords[1], 0.0]) }
            let minZ = (face .map { larger.vertices[$0].coords[2] }).min()!

            // face need to be above point for positive intersection

            if minZ <= maxZ {
                continue
            }
            
            let facebox = boundingBox( points:points )

            if boundingBoxInside( smaller:pointbox, larger:facebox ) {
                candidateFaces.append(fidx)
            }

        }

        if candidateFaces.count == 0 {
            return false 
        }

        // intersection if point x,y in face x,y

        for fidx in candidateFaces {

            let face = larger.faces[fidx]

            let verts = face .map { larger.vertices[$0] }
            let d1 = verts[1].sub(verts[0])
            let d2 = verts[2].sub(verts[0])

            let det = d1.coords[0]*d2.coords[1] - d1.coords[1]*d2.coords[0]

            if abs(det) < 0.000001 { continue }

            let dx = coord[0] - verts[0].coords[0]
            let dy = coord[1] - verts[0].coords[1]

            let r = (dx*d2.coords[1] - dy*d2.coords[0]) / det
            let s = (d1.coords[0]*dy - d1.coords[1]*dx) / det 

            if 0.0 <= r && r <= 1.0 && 0.0 <= s && s <= 1.0 {
                return true 
            }
        }

        return false 


    }

public func processNonMembraneTri( _ VERTICES:[Vector], _ NORMALS:[Vector],  _ FACES:[[Int]], _ opts:[String:Any]) 
                -> [(vertices:[Vector],normals:[Vector],faces:[[Int]],surfacetype:SurfaceType)]   {

    let COMPONENTS = surfaceComponentsDFS( faces:FACES, numvertices:VERTICES.count )

    var reentrantCOMPONENTS = [(vertices:[Vector],normals:[Vector],faces:[[Int]],surfacetype:SurfaceType)]()

    var probecenteredCOMPONENTS = [(vertices:[Vector],normals:[Vector],faces:[[Int]],surfacetype:SurfaceType)]()

    for component in COMPONENTS {

        if component.count < opts["minvertices"]! as! Int {
            continue
        }

        let vertices = component .map { VERTICES[$0] }
        let normals = component .map { NORMALS[$0] }

        // assume no boundary in non-membrane context

        let surftype = identifySurfaceType(vertices:vertices, normals:normals, boundaryvertices:[])

        if surftype == SurfaceType.undeterminedOpen {
            print("\nprocessNonMembraneTri : warning : have component with undetermined surface type, skipping")
            continue
        }

        if surftype == SurfaceType.reentrantClosed && opts["keepreentrant"]! as! Bool {

            let faces = findComponentFaces(component:component, FACES:FACES, numvertices:VERTICES.count)

            reentrantCOMPONENTS.append((vertices:vertices,normals:normals,faces:faces,surfacetype:surftype))
        }
        else if (opts["keepprobecentered"]! as! Bool) || !(opts["skipcavities"]! as! Bool) {

            let faces = findComponentFaces(component:component, FACES:FACES, numvertices:VERTICES.count)

            probecenteredCOMPONENTS.append((vertices:vertices,normals:normals,faces:faces,surfacetype:surftype))
        }

    }

    // sort in decreasing size 

    reentrantCOMPONENTS = reentrantCOMPONENTS .sorted { $0.vertices.count > $1.vertices.count }
    probecenteredCOMPONENTS = probecenteredCOMPONENTS .sorted { $0.vertices.count > $1.vertices.count }

    var keepCOMPONENTS = [(vertices:[Vector],normals:[Vector],faces:[[Int]],surfacetype:SurfaceType)]()

    if opts["keepreentrant"]! as! Bool {
        if reentrantCOMPONENTS.count > 0 {
            if opts["onlylargest"]! as! Bool {
                keepCOMPONENTS.append(reentrantCOMPONENTS[0])
            }
            else {
                keepCOMPONENTS += reentrantCOMPONENTS
            }
        }
    }

    if opts["keepprobecentered"]! as! Bool {
        if probecenteredCOMPONENTS.count > 0 {
            if opts["onlylargest"]! as! Bool {
                keepCOMPONENTS.append(probecenteredCOMPONENTS[0])
            }
            else {
                keepCOMPONENTS += probecenteredCOMPONENTS
            }
        }

        // since cavities are a subset of probe-centerd, just return

        return keepCOMPONENTS
    }


    var cavities = [Int]()

    // identify cavities

    if !(opts["skipcavities"]! as! Bool) && reentrantCOMPONENTS.count > 0 {

        // determine which reentrants enclose probe-centered, those are then assumed to be cavities

        for reent in reentrantCOMPONENTS {

            for (pidx,probecent) in probecenteredCOMPONENTS.enumerated() {
                if cavities.contains(pidx) || probecent.vertices.count > reent.vertices.count {
                    continue
                }

                if isInterior( smaller:probecent, larger:reent ) {
                    cavities.append(pidx)
                    probecenteredCOMPONENTS[pidx].surfacetype = SurfaceType.reentrantCavity
                }

            }
        }

        for pidx in cavities {
            keepCOMPONENTS.append(probecenteredCOMPONENTS[pidx])
        }

    }

    return keepCOMPONENTS

}









