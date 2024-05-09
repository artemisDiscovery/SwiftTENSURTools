
import Foundation
import MathTools

import Dispatch

let computeQueue = DispatchQueue( label:"compute", attributes: .concurrent )
let blocksQueue = DispatchQueue( label:"blocks" )


class AtomCircle {

    var center:Vector
    var radius:Double
    var exposure = [exposedArc]()
    var axis:AXES
    var atom:Int
    var removed:Bool
    var symmetryOp:String?

    init(_ atom:Int, _ center:Vector, _ radius:Double,  _ axis:AXES ) {
        self.atom = atom 
        self.center = center 
        self.radius = radius
        self.axis = axis
    
        removed = false

        symmetryOp = nil 
    }

    func setRemoved(_ state:Bool) {
        self.removed = state
    }


    func updateExposure( _ pstart:Vector, _ pend:Vector, _ othercircle: AtomCircle) {


        var newarc = exposedArc(pstart, pend, othercircle, othercircle, self )
        
        //print("updating exposure for atom \(atom)")
        //print("current count = \(exposure.count), newarc =  ")
        //print("\(newarc.str())")

        if exposure.count == 0 && !removed {
            exposure.append(newarc)
            return 
        }
  
        var newexposure = [exposedArc]()

        for iarc in 0..<exposure.count {
            newexposure += intersectArcs( exposure[iarc], newarc )
        }

        exposure = newexposure
        //Swift.print("new exposure atom \(self.atom) : ")
        //for arc in exposure {
        //    print(arc.str())
        //}

        if exposure.count == 0 {
            self.removed = true
            //Swift.print("in updateexposure, remove center at atom \(self.atom)")
        }


    }

    func intersectArcs(_ A: exposedArc, _ B: exposedArc ) -> [exposedArc] {
        
        let AstartInB = B.contains(A.ustart)
        let AendInB   = B.contains(A.uend)

        var returnArcs = [exposedArc]() 

        

        if AstartInB {
            if AendInB {
                // check for Astart, Aend in clockwise order in B
                if !B.clockwiseOrder(A.ustart, A.uend) {
                    var arc0 = exposedArc(B.pstart, A.pend, B.atomcircleStart, A.atomcircleEnd, self)
                    //arc0.setEndCircle(A.atomcircleEnd)
                    var arc1 = exposedArc(A.pstart, B.pend, A.atomcircleStart, B.atomcircleEnd, self)
                    //arc1.setEndCircle(B.atomcircleEnd)

                    returnArcs.append( arc0 )
                    returnArcs.append( arc1 )
                }
                else {
                    returnArcs.append( A )
                    
                }
                return returnArcs
            }
            else {
                // Astart in B and Aend NOT in B
                var arc0 = exposedArc(A.pstart, B.pend, A.atomcircleStart, B.atomcircleEnd, self)
                //arc0.setEndCircle(B.atomcircleEnd)
                returnArcs.append( arc0 )
                
                return returnArcs
            }
        }
        else if AendInB {
            // Astart NOT in B, Aend in B
            var arc0 = exposedArc(B.pstart, A.pend, B.atomcircleStart, A.atomcircleEnd, self)
            //arc0.setEndCircle(A.atomcircleEnd)
            returnArcs.append( arc0 )
            
            return returnArcs 
        }

        let BstartInA = A.contains(B.ustart)
        let BendInA   = A.contains(B.uend)

        if BstartInA {
            if BendInA {
                // check clockwise order Bstart, Bend
                if !A.clockwiseOrder( B.ustart, B.uend ) {
                    var arc0 = exposedArc(A.pstart, B.pend, A.atomcircleStart, B.atomcircleEnd, self)
                    //arc0.setEndCircle(B.atomcircleEnd)
                    var arc1 = exposedArc(B.pstart, A.pend, B.atomcircleStart, A.atomcircleEnd, self)
                    //arc1.setEndCircle(A.atomcircleEnd)
                    returnArcs.append( arc0 )
                    returnArcs.append( arc1 )

                }
                else {
                    returnArcs.append( B )
                }
                return returnArcs
            }
            else {
                // Bstart in A and Bend NOT in A
                var arc0 = exposedArc(B.pstart, A.pend, B.atomcircleStart, A.atomcircleEnd, self)
                //arc0.setEndCircle(A.atomcircleEnd)
                returnArcs.append( arc0 )
                return returnArcs
            }
        }
        else if BendInA {
            // Bstart NOT in A, Bend in A
            var arc0 = exposedArc(A.pstart, B.pend, A.atomcircleStart, B.atomcircleEnd, self)
            //arc0.setEndCircle(B.atomcircleEnd)
            returnArcs.append( arc0 )
            return returnArcs 
        }

        // no intersection

        return returnArcs 


    }

    
    public func str() -> String {
        var returnString = ""

        returnString += "center : \(center.coords[0]) , \(center.coords[1]), \(center.coords[2])\n"
        returnString += "radius : \(radius)\n"
        returnString += "axis : \(axis.rawValue)\n"
        returnString += "atom : \(atom)\n"
        returnString += "removed : \(removed)\n"
        returnString += "# exposed arcs : \(exposure.count)\n"

        for arc in exposure {
            returnString += arc.str()
        }

        return returnString
    }


    
}

public enum ContourError: Error {
    case missingCircleError
    case missingArcError
    
    
}


struct Contour {

    var axis:AXES
    var arcsInOrder = [exposedArc]() 
    var singletons = [AtomCircle]()
    var clockwise:Bool

    public init(_ circles: [AtomCircle]) throws  {
        // assume all have the same axis
        axis = circles[0].axis 

        var initialArc:exposedArc? = nil 
        var initialCircle:AtomCircle? = nil 

        for cidx in 0..<circles.count {
            var circle = circles[cidx]
            if circle.removed {
                continue 
            }

            if circle.exposure.count == 0 {
                singletons.append(circle)
                circle.setRemoved(true)
                continue
            }

            for arc in circle.exposure {
                if arc.removed {
                    continue
                }
                initialArc = arc 
                initialCircle = circle
                break
            }

            if initialArc != nil {
                break
            }

        }

        if initialArc == nil {
            print("contour: no initial arc")
            throw ContourError.missingArcError
        }


        var currentArc = initialArc!
        var currentCircle = initialCircle!

        //print("initial circle: ")
        //print(initialCircle!.str())
        //print("initial arc: ")
        //print(initialArc!.str())

        arcsInOrder.append(initialArc!)
        
        var closed = false
        var accumAngle = 0.0
        var prevdisp:Vector? = nil

        //while !closed {
        for _ in 0..<10000 {
            //print("contour: current arc :")
            //print(currentArc.str())
            let nextCircle = currentArc.atomcircleEnd
            
            if nextCircle.removed {
                print("contour: circle removed")
                //print("current circle :")
                //print(currentArc.parentcircle.str())
                //print("current arc :")
                //print(currentArc.str())
                //print("missing circle :")
                //print(nextCircle.str())
                throw ContourError.missingCircleError
            }

            var bestArc:exposedArc? = nil
            var bestDist:Double = Double.infinity

            //print("contour: next circle :")
            //print(nextCircle.str())

            for arc in nextCircle.exposure {
                //print("contour: test arc :")
                //print(arc.str())
                if arc.removed {
                    continue
                }
                let d = arc.pstart.sub(currentArc.pend).length()

                if d < bestDist {
                    bestArc = arc 
                    bestDist = d
                }

            }

            if bestArc == nil {
                print("contour: no best next arc found")
                throw ContourError.missingArcError
            }

            arcsInOrder.append(bestArc!)

            if bestArc!.pend == initialArc!.pstart {
                closed = true
                break
            }


            let disp = bestArc!.pstart.sub(currentArc.pstart).unit()
            if prevdisp != nil {
                var ang = acos(disp!.dot(prevdisp!))
                // define positive angle as clockwise
                // cross product expensive, better way?
                if disp!.cross(prevdisp!).dot(axisREF[axis.rawValue]) < 0.0 {
                    ang = -ang
                }

                accumAngle += ang
            }

            currentArc = bestArc!
            prevdisp = disp 
        }

        // assume closed here


        clockwise = accumAngle > 0.0

        for aidx in 0..<self.arcsInOrder.count {
            arcsInOrder[aidx].setRemoved(true)
        }

        
    }

    public func str(verbose:Bool=false) -> String {

        var report = ""

        report += "# singletons = \(singletons.count)\n"

        if verbose {
            for circle in singletons {
                report += circle.str()
            }
        }
        else {
            for circle in singletons {
                report += "\tatom \(circle.atom)\n"
            }
        }

        report += "\n# arcs = \(arcsInOrder.count)\n"

        if verbose {
            for arc in arcsInOrder {
                report += arc.str()
            }
        }
        else {
            for arc in arcsInOrder {
                report += "\tstart atom : \(arc.atomcircleStart.atom), end atom : \(arc.atomcircleEnd.atom)\n"
            }
        }

        report += "\nclockwise : \(clockwise)"

        return report 

    }


}

// ref axis is just unit vector along X, Y or Z

let axisREF = [Vector([1.0,0.0,0.0]), Vector([0.0,1.0,0.0]), Vector([0.0,0.0,1.0])]

// 'right' and 'up' are references for planes perpendicular to X Y OR Z
// X -> Y,Z, Y -> Z,X , Z -> X,Y

let axisRIGHT = [1, 2, 0]
let axisUP    = [2, 0, 1]


struct exposedArc {

    // pstart, pend points to define

    var acute:Bool 
    var dotSE:Double
    var clockwiseMeasure:Double
    var pstart:Vector
    var pend:Vector

    var ustart:Vector
    var uend:Vector


    var parentcircle:AtomCircle
    

    // intersections with 
    var atomcircleStart:AtomCircle
    var atomcircleEnd:AtomCircle

    var removed:Bool

    
    init( _ start:Vector, _ end:Vector, _ startcircle: AtomCircle, _ endcircle: AtomCircle, _ parentcircle:AtomCircle ) {
        
            pstart = start
            pend = end 
            atomcircleStart = startcircle
            atomcircleEnd = endcircle

            self.parentcircle = parentcircle

            //self.axis = parentcircle.axis
            //self.center = parentcircle.center
            


            ustart = start.sub(parentcircle.center)
            ustart = ustart.scale(1.0/ustart.length())

            uend = end.sub(parentcircle.center)
            uend = uend.scale(1.0/uend.length())

            let ref = axisREF[parentcircle.axis.rawValue]

            acute = ref.dot(uend.cross(ustart)) > 0.0

            dotSE = ustart.dot(uend)

            // I can't call the function ??

            // clockwiseMeasure = clockwise(uend)

            let crossdot = uend.cross(ustart).dot(ref)

            if crossdot >= 0 {
                clockwiseMeasure =  1.0 - uend.dot(ustart)
            }
            else {
                clockwiseMeasure =  3.0 + uend.dot(ustart)
            }

                removed = false 
 

    }

    mutating func setEndCircle( _ end: AtomCircle) {
        self.atomcircleEnd = end
    }

    mutating func setRemoved(_ state:Bool) {
        self.removed = state 
    }

    // return inside/outside for two vectors

    func contains( _ u:Vector) -> Bool {
        

        let dotS = ustart.dot(u)
        let dotE = uend.dot(u)

        // cosine monotonic - if max is greater than dotSE, then vector is inside smaller angle
        // BUT, this fast test ONLY WORKS for 'smaller angle' < 120 degrees; otherwise vector at (360 - theta)/2 can violate test

        // 
        if dotSE > -0.5 {

            if min(dotS,dotE) > dotSE {
                // inside smaller angle - note that 'acute' is not properly used here, I meant '< 180 deg'
                if acute {
                    return true
                }
                else {
                    return false
                }
            }
            else {
                // outside smaller angle
                if acute {
                    return false 
                }
                else {
                    return true
                }
            }
        }
        
        // need to do slow test

        return clockwise(u) < clockwiseMeasure

        
    }

    func clockwise( _ u:Vector  ) -> Double {
        // return positive measure of clockwise angle w.r.t start of arc

        let ref = axisREF[parentcircle.axis.rawValue]

        let crossdot = u.cross(ustart).dot(ref)

        if crossdot >= 0 {
            return 1.0 - u.dot(ustart)
        }
        else {
            return 3.0 + u.dot(ustart)
        }

    }

    func clockwiseOrder( _ uA:Vector, _ uB:Vector ) -> Bool {

        let measureA = clockwise( uA )
        let measureB = clockwise( uB )

        return measureA < measureB 
    }

    public func str() -> String {
        var returnString = ""

        returnString += "\tarc pstart = \(pstart.coords[0]), \(pstart.coords[1]), \(pstart.coords[2])\n"
        returnString += "\tarc pend   = \(pend.coords[0]), \(pend.coords[1]), \(pend.coords[2])\n"
        returnString += "\tarc ustart = \(ustart.coords[0]), \(ustart.coords[1]), \(ustart.coords[2])\n"
        returnString += "\tarc uend   = \(uend.coords[0]), \(uend.coords[1]), \(uend.coords[2])\n"
        returnString += "\tparent circle atom = \(parentcircle.atom)\n"
        returnString += "\tstart intersection atom = \(atomcircleStart.atom)\n"
        returnString += "\tend intersection atom = \(atomcircleEnd.atom)\n"
        returnString += "\tdotSE : \(dotSE)\n"
        returnString += "\tacute : \(acute)\n"
        returnString += "\tremoved : \(removed)\n"
        


        return returnString
    }
    
}

// find intersection between circles, update exposed arcs and/or remove buried circle

func intersectAtomCircles(_ circleA: AtomCircle, _ circleB: AtomCircle) {

    if circleA.removed || circleB.removed {
        return
    }

    var uAB = circleB.center.sub(circleA.center)
    let dAB = uAB.length()
    let rA = circleA.radius
    let rB = circleB.radius
    // no intersection 

    if dAB > (rA + rB) {
        //print("no intersection")
        return
    }

    // B in A 
    if dAB + circleB.radius < circleA.radius {
        circleB.setRemoved(true) 
        //print("remove circle \(circleB.atom), inside \(circleA.atom)")
        return
    }

    // A in B 
    if dAB + circleA.radius < circleB.radius {
        circleA.setRemoved(true)
        //print("remove circle \(circleA.atom), inside \(circleB.atom)")
        return
    }

    let RIGHT = axisRIGHT[circleA.axis.rawValue]
    let UP = axisUP[circleA.axis.rawValue]
    uAB = uAB.scale(1.0/dAB)
    var vAB = Vector([0.0,0.0,0.0])
    vAB.coords[RIGHT] = -uAB.coords[UP]
    vAB.coords[UP] = uAB.coords[RIGHT]

    let m = (dAB*dAB + rA*rA - rB*rB)/(2.0 * dAB)

    // check for valid intersection

    let disc = rA*rA - m*m 

    if disc <= 0.0 {
        return 
    }

    let h = sqrt(disc)

    let pdown = circleA.center.add(uAB.scale(m)).add(vAB.scale(-h))
    let pup = circleA.center.add(uAB.scale(m)).add(vAB.scale(h))

    // clockwise ; for A, pdown to pup, for B, pup to pdown

    circleA.updateExposure(pdown, pup, circleB)
    circleB.updateExposure(pup, pdown, circleA)


}



func circlesForAtoms( atompos:Matrix<Double>, radii:[Double], proberad:Double,
       axis:AXES, minCoord:Double, delta:Double, limits:[Int], thread:Int ) 
        -> ([(AtomCircle,Int)],Int)   {
        
    var data = [(AtomCircle,Int)]()
    let storage = atompos.getStorage()
    let RIGHT = axisRIGHT[axis.rawValue]
    let UP = axisUP[axis.rawValue]

    for aidx in limits[0]..<limits[1] {
        
        let catom = Vector((0..<3) .map { storage[3*aidx + $0]})
        let caxis = catom.coords[axis.rawValue]

        let augrad = radii[aidx] + proberad

        // bottom layer is below bottom of augmented sphere, so start one higher
        let layerlo = Int(((caxis - augrad) - minCoord)/delta) + 1
        // upper layer contains top of augmented sphere, include it in range
        let layerhi = Int(((caxis + augrad) - minCoord)/delta) + 1

        for layer in layerlo..<layerhi {
            let layerc = minCoord + Double(layer)*delta 
            let h = abs(caxis - layerc) 
            // sanity test
            if h > augrad {
                continue
            }
            let crad = sqrt(augrad*augrad - h*h)

            var center = Vector([0.0,0.0,0.0])
            center.coords[axis.rawValue] = layerc
            center.coords[RIGHT] = catom.coords[RIGHT]
            center.coords[UP] = catom.coords[UP]
            
            data.append((AtomCircle(aidx, center, crad, axis ),layer))
        }
    }

    return (data,thread)
        
}

func addBLOCK( _ BLOCKS: inout [[(AtomCircle,Int)]?], _ data:([(AtomCircle,Int)],Int) ) {
    BLOCKS[data.1] = data.0
}

func atomCirclesForLayers( atompos:Matrix<Double>, radii:[Double], 
    proberad:Double, minaxiscoord:Double, layerdelta:Double, axis:AXES, numthreads:Int=1 ) -> LAYERS {

    
    let shape = atompos.getShape()
    let numatoms = shape[0]
    let storage = atompos.getStorage()
    let axiscoords = (0..<numatoms) .map { storage[3*$0 + axis.rawValue] }


    var LIMITS = [[Int]]()

    let atomsPerChunk = Int(shape[0]/numthreads)
    let nchunks = Int(ceil(Double(shape[0])/Double(atomsPerChunk)))

    for ichunk in 0..<nchunks {
        let lo = ichunk*atomsPerChunk
        var hi = lo + atomsPerChunk
        if ichunk == nchunks - 1 {
            hi = numatoms
        }

        LIMITS.append([lo,hi])
    }  

    // each block will hold list of atom circles coupled with layer index

    var BLOCKS:[[(AtomCircle,Int)]?] = Array(repeating:nil, count:nchunks)

    let group = DispatchGroup() 

    for ichunk in 0..<nchunks {

        group.enter() 

        computeQueue.async {
            let data = circlesForAtoms( atompos:atompos, radii:radii, proberad:proberad,
                axis:axis, minCoord:minaxiscoord, delta:layerdelta, limits:LIMITS[ichunk], 
                thread:ichunk ) 

            blocksQueue.sync {
                addBLOCK( &BLOCKS, data )
            }

            group.leave()

        }

        

    }

    group.wait()

    var circles = [AtomCircle]()
    var circleLayers = [Int]() 

    for block in BLOCKS {
        circles += block! .map { $0.0 }
        circleLayers += block! .map { $0.1 }
    }

    let circleLAYERS = LAYERS(circles, circleLayers, layerdelta, axis, minaxiscoord )

    return circleLAYERS
}

func intersectCirclesInLayerRange( _ circleLayers:LAYERS, _ limits:[Int]) -> ([[AtomCircle]],[Int])? {

    let layerBits = circleLayers.layerBits

    var allcircles = circleLayers.objects as! [AtomCircle]

    var returnCircles = [[AtomCircle]]()
    var returnLayerIndices = [Int]()

    for lidx in limits[0]..<limits[1] {
        if layerBits[lidx] == nil {
            continue
        }

        var layercircles = layerBits[lidx]!.indices() .map { allcircles[$0] }

        if layercircles.count == 1 {
            returnCircles.append(layercircles)
            returnLayerIndices.append(lidx)
            continue
        }

        var centers = [Double]()
        var radii = [Double]()

        for circle in layercircles {
            centers += circle.center.coords
            radii.append(circle.radius)
        }

        let circleCenters = Matrix<Double>([layercircles.count,3], content:centers )
        let circleRadii = Matrix<Double>([layercircles.count], content:radii )

        var pairs:[[Int]]?

        do {
            var dists = try cdist(circleCenters,circleCenters)
            try dists.setdiagonal(1.0e12)
            let sumRadii = try circleRadii.addTranspose(circleRadii)

            let mask = try Mask.compare(dists, sumRadii) { $0 < $1 }
            pairs = mask.nonzero()
        }
        catch {
            print("exception in intersectCircles for layer \(lidx)")
            return nil 
        }


        for p in pairs! {
            if p[0] < p[1] {
                intersectAtomCircles( layercircles[p[0]], layercircles[p[1]] )

            }
        }

        returnCircles.append(layercircles)
        returnLayerIndices.append(lidx)


    }

    return (returnCircles,returnLayerIndices)

}

func intersectingCirclesForLayers( _ circleLayers:LAYERS, numthreads:Int=1 ) {

    // set limits using fraction of total number of circles

    var LIMITS = [[Int]]()

    let atomcircles = circleLayers.objects as! [AtomCircle]
    let layerBits = circleLayers.layerBits

    let circlesPerChunk = Int(atomcircles.count/numthreads)

    // lo and hi specify layer ranges

    var currentLo = 0
    var currentCircleCount = 0
     
    while currentLo < layerBits.count  {

        
        if layerBits[currentLo] != nil {
            currentCircleCount = layerBits[currentLo]!.indices().count
        }


        if currentLo == layerBits.count - 1 {
            break 
        }

        var lastHi = -1 

        for currentHi in (currentLo+1)..<layerBits.count {
            if layerBits[currentHi] != nil {
                currentCircleCount += layerBits[currentHi]!.indices().count
            }
            if currentCircleCount >= circlesPerChunk {
                LIMITS.append([currentLo,currentHi+1])
                currentCircleCount = 0 
                lastHi = currentHi
                break
            }
        }

        if lastHi < 0 {
            break
        }

        currentLo = lastHi + 1

    }

    // add last if needed 

    if currentCircleCount > 0 {
        LIMITS.append([currentLo,layerBits.count])
    }

    var cumcount = 0 

    for lim in LIMITS {
        print("limit \(lim)")
        var ccount = 0 
        for lidx in lim[0]..<lim[1] {
            if layerBits[lidx] == nil {
                continue 
            }
            ccount += layerBits[lidx]!.indices().count
        }
        print("\tcount for range = \(ccount)")
        cumcount += ccount
    }

    print("\ncumulative count = \(cumcount), compare to \(atomcircles.count) circles expected")

    var BLOCKS = [([[AtomCircle]],[Int])]()

    let group = DispatchGroup() 

    for ichunk in 0..<LIMITS.count {

        group.enter() 

        computeQueue.async  {
            let data = intersectCirclesInLayerRange(circleLayers,LIMITS[ichunk]) 

            blocksQueue.sync {
                BLOCKS.append(data!)
            }

            group.leave()

        }

        

    }

    group.wait()

    BLOCKS = BLOCKS.sorted { $0.1[0] < $1.1[0] }

    var sortedLayers = [Int]()
    var sortedCircles = [[AtomCircle]]()

    for block in BLOCKS {
        sortedLayers += block.1
        sortedCircles += block.0
    }

    for lidx in 0..<sortedLayers.count {
        print("sorted layer \(sortedLayers[lidx]), \(sortedCircles[lidx].count) circles")
        
    }

    var contoursForLevels = [[Contour]]()

    for lidx in 0..<sortedLayers.count {
        var contours = [Contour]()
        print("layer \(lidx) contours ...")
        for j in 0..<2 {
            var cont:Contour?
            do {
                cont = try Contour(sortedCircles[lidx])
            }
            catch {
                print("exception in Contour")
                break
            }
            
            if cont == nil {
                print("nil contour, continue")
                continue
            }
            if cont!.singletons.count == 0 && cont!.arcsInOrder.count == 0 {
                break
            }
            contours.append(cont!)

        }
        contoursForLevels.append(contours)
        print("layer \(lidx) \(contours.count) contours ...")
        //for cont in contours {
        //    print(cont.str())
        //}
    }

    


}