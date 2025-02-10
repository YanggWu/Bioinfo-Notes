hljs.registerLanguage('shell', function(hljs) {
  var KEYWORDS = {
    keyword: 'if then else fi for while in do done exit return break continue case esac',
    literal: 'true false',
    built_in: 'echo read cd pwd mkdir git gemma bwa hisat2 samtools bcftools'  // 添加外部命令到内置命令列表
  };
  return {
    keywords: KEYWORDS,
    contains: [
      hljs.HASH_COMMENT_MODE,
      hljs.BACKSLASH_ESCAPE
    ]
  };
});

hljs.initHighlightingOnLoad();