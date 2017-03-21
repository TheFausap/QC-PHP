<?php

require_once __DIR__ . '/vendor/autoload.php';
require_once 'complex.php';

use MathPHP\LinearAlgebra\Matrix;
use MathPHP\LinearAlgebra\Vector;

function tensorProduct(Vector $v1, Vector $v2)
{
	$lv1 = $v1->getN();
	$lv2 = $v2->getN();
	$k = 0;
	$ra = array_fill(0,$lv1*$lv2-1,0.0);

	for ($i = 0; $i < $lv1; $i++) {
		for ($j = 0; $j < $lv2; $j++) {
			$ra[$k] = $v1[$i]*$v2[$j];
			$k += 1;
		}
	}			
	$r = new Vector($ra);
	return $r;
}

$k1 = new Vector([0.0, 1.0]);
$k0 = new Vector([1.0, 0.0]);
$k00 = tensorProduct($k0,$k0);

$I = complex(0.0,1.0);
$MI = complex(0.0,-1.0);

$invsqr2 = 1.0/sqrt(2.0);

$H = new Matrix(
	[
		[1.0,  1.0],
		[1.0, -1.0]
	]
);

$Y = new Matrix(
	[
		[0.0,  $I],
		[$MI, 0.0]
	]
);

$H = $H->scalarMultiply($invsqr2);
$HH = $H->kroneckerProduct($H);

$prod  = $H->vectorMultiply($k0);
//$prod1 = $Y->vectorMultiply($k0);

function printq(Vector $v, int $width)
{
	$len = $v->getN();

	$noexpand = false;

	if ($width == 0) {
		$noexpand = true;
	}

	if ($noexpand) {
		for ($i = 0; $i < $len; $i++) {
			echo $v->get($i) . "|" . decbin($i) . "> ";
		}
	}	

	echo "\n";
}

printq($prod,0);
echo $k00 . " done \n";

echo $HH . " done \n";

//echo "Vector: " . $prod1 . "\n";

?>
